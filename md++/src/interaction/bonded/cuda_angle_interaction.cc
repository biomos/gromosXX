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

  m_initialized = true;

  if (!quiet)
    os << "CUDA ANGLE INTERACTION\n"
       << "\tterms: " << m_num_angles << "\n"
       << "END\n";
  return 0;
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
  cudaMemset(m_angle_energy.data(), 0, sizeof(double) * num_energy_groups);
  cudaMemset(m_virial.data(), 0, sizeof(double) * 9);

  const math::CuVArray::View pos =
      sim.cuda().configuration_view(conf, gpu::MIRROR_POS).current().pos;

  static gpu::cuvector<FPL3_TYPE> force;
  const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
  if (force.size() < num_atoms) force.resize(num_atoms);
  cudaMemset(force.data(), 0, sizeof(FPL3_TYPE) * num_atoms);

  gpu::launch_angle(
      pos, m_angle_i.data(), m_angle_j.data(), m_angle_k.data(), m_angle_type.data(),
      m_K.data(), m_cos0.data(), m_atom_energy_group.data(), m_num_angles,
      conf.boundary_type, conf.current().box,
      force.data(), m_angle_energy.data(), m_virial.data());

  cudaDeviceSynchronize();

  for (unsigned i = 0; i < num_atoms; ++i) {
    conf.current().force(i) += math::Vec(force[i].x, force[i].y, force[i].z);
  }

  for (unsigned g = 0; g < num_energy_groups; ++g) {
    conf.current().energies.angle_energy[g] += m_angle_energy[g];
  }

  for (unsigned b = 0; b < 3; ++b) {
    for (unsigned a = 0; a < 3; ++a) {
      conf.current().virial_tensor(b, a) += m_virial[b * 3 + a];
    }
  }

  m_timer.stop();
  return 0;
}
