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
 * @file cuda_position_restraint_interaction.cc
 * GPU-native position-restraint interaction. See
 * cuda_position_restraint_interaction.h. Plain C++, no kernel syntax;
 * the real __global__ kernel lives in
 * gpu/cuda/interaction/special/position_restraint_kernels.cu, compiled
 * into grocuda.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"
#include "../../interaction/interaction_types.h"

#include "../../gpu/cuda/manager/cuda_manager.h"
#include "gpu/cuda/interaction/special/position_restraint_kernels.h"
#include "gpu/cuda/memory/energy_accumulate_kernels.h"

#include "cuda_position_restraint_interaction.h"

#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

int interaction::CUDA_Position_Restraint_Interaction::init(
    topology::Topology & topo,
    configuration::Configuration & conf,
    simulation::Simulation & sim,
    std::ostream & os,
    bool quiet) {

  const std::vector<topology::position_restraint_struct> & restraints =
      topo.position_restraints();
  m_num_restraints = static_cast<unsigned>(restraints.size());

  m_seq.resize(m_num_restraints);
  m_ref.resize(m_num_restraints);
  m_inv_bf_scale.resize(m_num_restraints);
  m_atom_energy_group.resize(m_num_restraints);

  const bool use_bfactor = sim.param().posrest.posrest == simulation::posrest_bfactor;
  const math::VArray & ref = conf.special().reference_positions;
  const math::SArray & bfactor = conf.special().bfactors;

  for (unsigned t = 0; t < m_num_restraints; ++t) {
    const unsigned a = static_cast<unsigned>(restraints[t].seq);
    m_seq[t] = a;
    m_ref[t] = FPL3_TYPE{static_cast<FPL_TYPE>(ref(a)(0)),
                          static_cast<FPL_TYPE>(ref(a)(1)),
                          static_cast<FPL_TYPE>(ref(a)(2))};
    m_inv_bf_scale[t] = use_bfactor
        ? static_cast<FPL_TYPE>(1.0 / bfactor(a))
        : static_cast<FPL_TYPE>(1.0);
    m_atom_energy_group[t] = topo.atom_energy_group()[a];
  }

  const unsigned num_energy_groups = static_cast<unsigned>(topo.energy_groups().size());
  m_posrest_energy.resize(num_energy_groups);

  if (m_stream == 0) cudaStreamCreate(&m_stream);

  m_initialized = true;

  if (!quiet)
    os << "CUDA POSITION RESTRAINT INTERACTION\n"
       << "\tterms: " << m_num_restraints << "\n"
       << "END\n";
  return 0;
}

interaction::CUDA_Position_Restraint_Interaction::~CUDA_Position_Restraint_Interaction() {
  if (m_stream) cudaStreamDestroy(m_stream);
}

int interaction::CUDA_Position_Restraint_Interaction::calculate_interactions(
    topology::Topology & topo,
    configuration::Configuration & conf,
    simulation::Simulation & sim) {

  m_timer.start(sim);

  if (!m_initialized) {
    m_timer.stop();
    return 1;
  }

  const unsigned num_energy_groups = static_cast<unsigned>(m_posrest_energy.size());
  cudaMemsetAsync(m_posrest_energy.data(), 0, sizeof(double) * num_energy_groups, m_stream);

  // Force written directly into the GPU-resident mirror (zeroed once per
  // step by Forcefield::calculate_interactions()'s sim.cuda().
  // zero_mirror_force(), not here) -- no private scratch buffer, no
  // sync, no host readback. Only requesting MIRROR_POS as a read field
  // (never MIRROR_FORCE) is what keeps this from triggering a coarse
  // resync that would clobber force other GPU-native terms may have
  // already accumulated into the mirror this step.
  gpu::Configuration::View view = sim.cuda().configuration_view(conf, gpu::MIRROR_POS, m_stream);

  if (!m_energy_registered) {
    sim.cuda().ensure_energy_groups(conf, num_energy_groups);
    m_energy_registered = true;
  }

  gpu::launch_position_restraint(
      view.current().pos, m_seq.data(), m_ref.data(), m_inv_bf_scale.data(),
      static_cast<FPL_TYPE>(sim.param().posrest.force_constant),
      m_atom_energy_group.data(), m_num_restraints,
      conf.boundary_type, conf.current().box,
      view.current().force.data(), m_posrest_energy.data(), m_stream);

  sim.cuda().mark_gpu_dirty(conf, gpu::MIRROR_FORCE, m_stream);

  // Same on-device merge for energy -- see energy_accumulate_kernels.h
  // and cuda_angle_interaction.cc's identical comment. No host sync at
  // all in this class anymore.
  const gpu::EnergyMirrorPtrs eptrs = sim.cuda().energy_mirror_ptrs(conf);
  gpu::launch_accumulate_energy(eptrs.posrest, m_posrest_energy.data(), num_energy_groups, m_stream);
  sim.cuda().mark_gpu_dirty(conf, gpu::MIRROR_ENERGY, m_stream);

  m_timer.stop();
  return 0;
}
