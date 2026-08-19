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
 * @file cuda_improper_dihedral_interaction.cc
 * GPU-native improper dihedral interaction. See
 * cuda_improper_dihedral_interaction.h. Plain C++, no kernel syntax; the
 * real __global__ kernel lives in
 * gpu/cuda/interaction/bonded/improper_dihedral_kernels.cu, compiled
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
#include "gpu/cuda/interaction/bonded/improper_dihedral_kernels.h"

#include "cuda_improper_dihedral_interaction.h"

#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE bonded

int interaction::CUDA_Improper_Dihedral_Interaction::init(
    topology::Topology & topo,
    configuration::Configuration & conf,
    simulation::Simulation & sim,
    std::ostream & os,
    bool quiet) {

  if (sim.param().perturbation.perturbation) {
    io::messages.add(
        "CUDA_Improper_Dihedral_Interaction does not support perturbation "
        "(out of scope, matching every other perturbation gate in the "
        "CUDA path).",
        "CUDA_Improper_Dihedral_Interaction", io::message::error);
    return 1;
  }
  if (sim.param().gamd.gamd) {
    io::messages.add(
        "CUDA_Improper_Dihedral_Interaction does not support GAMD.",
        "CUDA_Improper_Dihedral_Interaction", io::message::error);
    return 1;
  }

  const std::vector<topology::four_body_term_struct> & dihedrals =
      topo.solute().improper_dihedrals();
  m_num_dihedrals = static_cast<unsigned>(dihedrals.size());

  m_dihedral_i.resize(m_num_dihedrals);
  m_dihedral_j.resize(m_num_dihedrals);
  m_dihedral_k.resize(m_num_dihedrals);
  m_dihedral_l.resize(m_num_dihedrals);
  m_dihedral_type.resize(m_num_dihedrals);
  m_atom_energy_group.resize(m_num_dihedrals);

  for (unsigned t = 0; t < m_num_dihedrals; ++t) {
    m_dihedral_i[t] = static_cast<unsigned>(dihedrals[t].i);
    m_dihedral_j[t] = static_cast<unsigned>(dihedrals[t].j);
    m_dihedral_k[t] = static_cast<unsigned>(dihedrals[t].k);
    m_dihedral_l[t] = static_cast<unsigned>(dihedrals[t].l);
    m_dihedral_type[t] = static_cast<unsigned>(dihedrals[t].type);
    m_atom_energy_group[t] = topo.atom_energy_group()[dihedrals[t].i];
  }

  const std::vector<interaction::improper_dihedral_type_struct> & types =
      topo.impdihedral_types();
  m_K.resize(types.size());
  m_q0.resize(types.size());
  for (unsigned t = 0; t < types.size(); ++t) {
    m_K[t] = static_cast<FPL_TYPE>(types[t].K);
    m_q0[t] = static_cast<FPL_TYPE>(types[t].q0);
  }

  const unsigned num_energy_groups = static_cast<unsigned>(topo.energy_groups().size());
  m_improper_energy.resize(num_energy_groups);
  m_virial.resize(9);

  if (m_stream == 0) cudaStreamCreate(&m_stream);

  m_initialized = true;

  if (!quiet)
    os << "CUDA IMPROPER DIHEDRAL INTERACTION\n"
       << "\tterms: " << m_num_dihedrals << "\n"
       << "END\n";
  return 0;
}

interaction::CUDA_Improper_Dihedral_Interaction::~CUDA_Improper_Dihedral_Interaction() {
  if (m_stream) cudaStreamDestroy(m_stream);
}

int interaction::CUDA_Improper_Dihedral_Interaction::calculate_interactions(
    topology::Topology & topo,
    configuration::Configuration & conf,
    simulation::Simulation & sim) {

  m_timer.start(sim);

  if (!m_initialized) {
    m_timer.stop();
    return 1;
  }

  const unsigned num_energy_groups = static_cast<unsigned>(m_improper_energy.size());
  cudaMemsetAsync(m_improper_energy.data(), 0, sizeof(double) * num_energy_groups, m_stream);
  cudaMemsetAsync(m_virial.data(), 0, sizeof(double) * 9, m_stream);

  gpu::Configuration::View view = sim.cuda().configuration_view(conf, gpu::MIRROR_POS);

  gpu::launch_improper_dihedral(
      view.current().pos, m_dihedral_i.data(), m_dihedral_j.data(), m_dihedral_k.data(),
      m_dihedral_l.data(), m_dihedral_type.data(),
      m_K.data(), m_q0.data(), m_atom_energy_group.data(), m_num_dihedrals,
      conf.boundary_type, conf.current().box,
      view.current().force.data(), m_improper_energy.data(), m_virial.data(), m_stream);

  sim.cuda().mark_gpu_dirty(conf, gpu::MIRROR_FORCE);

  cudaStreamSynchronize(m_stream);

  for (unsigned g = 0; g < num_energy_groups; ++g) {
    conf.current().energies.improper_energy[g] += m_improper_energy[g];
  }

  for (unsigned b = 0; b < 3; ++b) {
    for (unsigned a = 0; a < 3; ++a) {
      conf.current().virial_tensor(b, a) += m_virial[b * 3 + a];
    }
  }

  m_timer.stop();
  return 0;
}
