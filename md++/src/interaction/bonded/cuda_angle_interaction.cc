/*
 * This file is part of GROMOS.
 *
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 *
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file cuda_angle_interaction.cc
 * GPU-native harmonic (cosine) bond-angle interaction. See
 * cuda_angle_interaction.h. Plain C++, no kernel syntax; the real
 * __global__ kernel lives in
 * gpu/cuda/interaction/bonded/angle_kernels.cu, compiled into grocuda.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"
#include "../../interaction/interaction_types.h"

#include "../../gpu/cuda/manager/cuda_manager.h"
#include "gpu/cuda/interaction/bonded/angle_kernels.h"
#include "gpu/cuda/memory/virial_accumulate_kernels.h"

#include "cuda_angle_interaction.h"

#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE bonded

int interaction::CUDA_Angle_Interaction::init(
    topology::Topology & topo,
    configuration::Configuration & conf,
    simulation::Simulation & sim,
    std::ostream & os,
    bool quiet) {

  if (sim.param().perturbation.perturbation) {
    io::messages.add(
        "CUDA_Angle_Interaction does not support perturbation "
        "(out of scope, matching every other perturbation gate in the "
        "CUDA path).",
        "CUDA_Angle_Interaction", io::message::error);
    return 1;
  }
  if (sim.param().gamd.gamd) {
    io::messages.add(
        "CUDA_Angle_Interaction does not support GAMD.",
        "CUDA_Angle_Interaction", io::message::error);
    return 1;
  }

  const std::vector<topology::three_body_term_struct> & angles = topo.solute().angles();
  m_num_angles = static_cast<unsigned>(angles.size());

  m_angle_i.resize(m_num_angles);
  m_angle_j.resize(m_num_angles);
  m_angle_k.resize(m_num_angles);
  m_angle_type.resize(m_num_angles);
  m_atom_energy_group.resize(m_num_angles);

  for (unsigned t = 0; t < m_num_angles; ++t) {
    m_angle_i[t] = static_cast<unsigned>(angles[t].i);
    m_angle_j[t] = static_cast<unsigned>(angles[t].j);
    m_angle_k[t] = static_cast<unsigned>(angles[t].k);
    m_angle_type[t] = static_cast<unsigned>(angles[t].type);
    m_atom_energy_group[t] = topo.atom_energy_group()[angles[t].i];
  }

  const std::vector<interaction::angle_type_struct> & angletypes = topo.angle_types_cosharm();
  m_K.resize(angletypes.size());
  m_cos0.resize(angletypes.size());
  for (unsigned t = 0; t < angletypes.size(); ++t) {
    m_K[t] = static_cast<FPL_TYPE>(angletypes[t].K);
    m_cos0[t] = static_cast<FPL_TYPE>(angletypes[t].cos0);
  }

  const unsigned num_energy_groups = static_cast<unsigned>(topo.energy_groups().size());
  m_angle_energy.resize(num_energy_groups);
  m_virial.resize(9);

  if (m_stream == 0) cudaStreamCreate(&m_stream);

  m_initialized = true;

  if (!quiet)
    os << "CUDA ANGLE INTERACTION\n"
       << "\tterms: " << m_num_angles << "\n"
       << "END\n";
  return 0;
}

interaction::CUDA_Angle_Interaction::~CUDA_Angle_Interaction() {
  if (m_stream) cudaStreamDestroy(m_stream);
}

int interaction::CUDA_Angle_Interaction::calculate_interactions(
    topology::Topology & topo,
    configuration::Configuration & conf,
    simulation::Simulation & sim) {

  m_timer.start(sim);

  if (!m_initialized) {
    m_timer.stop();
    return 1;
  }

  const unsigned num_energy_groups = static_cast<unsigned>(m_angle_energy.size());
  cudaMemsetAsync(m_angle_energy.data(), 0, sizeof(double) * num_energy_groups, m_stream);
  cudaMemsetAsync(m_virial.data(), 0, sizeof(double) * 9, m_stream);

  // Force written directly into the GPU-resident mirror (zeroed once
  // per step by Forcefield::calculate_interactions()'s sim.cuda().
  // zero_mirror_force(), not here) -- no private scratch buffer, no
  // sync, no host readback. Only requesting MIRROR_POS as a read field
  // (never MIRROR_FORCE) is what keeps this from triggering a coarse
  // resync that would clobber the force NonBonded/other bonded terms
  // may have already accumulated this step.
  gpu::Configuration::View view = sim.cuda().configuration_view(conf, gpu::MIRROR_POS, m_stream);

  gpu::launch_angle(
      view.current().pos, m_angle_i.data(), m_angle_j.data(), m_angle_k.data(), m_angle_type.data(),
      m_K.data(), m_cos0.data(), m_atom_energy_group.data(), m_num_angles,
      conf.boundary_type, conf.current().box,
      view.current().force.data(), m_angle_energy.data(), m_virial.data(), m_stream);

  sim.cuda().mark_gpu_dirty(conf, gpu::MIRROR_FORCE, m_stream);

  // Publish virial into the shared mirror's virial_tensor via atomicAdd
  // (GPU-resident, no CPU round trip) -- every bonded term/NonBonded/
  // active constraint algorithm contributes to the same global
  // accumulator each step, so unlike constraint_force this genuinely
  // needs atomicAdd, not a plain write. Zeroed once per step by
  // CudaManager::zero_mirror_force().
  gpu::launch_accumulate_virial9(
      reinterpret_cast<FPH_TYPE*>(view.current().virial_tensor), m_virial.data(), m_stream);
  sim.cuda().mark_gpu_dirty(conf, gpu::MIRROR_VIRIAL, m_stream);

  // Energy is a small, private, double-precision buffer (not part of
  // the mirror -- see cuda_angle_interaction.h) so still needs a sync
  // to read back, but only on this algorithm's own stream, not a
  // device-wide barrier that would stall NonBonded/other bonded terms
  // running concurrently.
  cudaStreamSynchronize(m_stream);

  for (unsigned g = 0; g < num_energy_groups; ++g) {
    conf.current().energies.angle_energy[g] += m_angle_energy[g];
  }

  m_timer.stop();
  return 0;
}
