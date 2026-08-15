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
 * @file cuda_dihedral_interaction.cc
 * GPU-native torsional dihedral interaction. See
 * cuda_dihedral_interaction.h. Plain C++, no kernel syntax; the real
 * __global__ kernel lives in
 * gpu/cuda/interaction/bonded/dihedral_kernels.cu, compiled into
 * grocuda.
 */

#include "../../stdheader.h"

#include <array>

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"
#include "../../interaction/interaction_types.h"

#include "../../gpu/cuda/manager/cuda_manager.h"
#include "gpu/cuda/interaction/bonded/dihedral_kernels.h"
#include "gpu/cuda/interaction/bonded/sparse_force_accumulate.h"

#include "cuda_dihedral_interaction.h"

#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE bonded

int interaction::CUDA_Dihedral_Interaction::init(
    topology::Topology & topo,
    configuration::Configuration & conf,
    simulation::Simulation & sim,
    std::ostream & os,
    bool quiet) {

  if (sim.param().perturbation.perturbation) {
    io::messages.add(
        "CUDA_Dihedral_Interaction does not support perturbation "
        "(out of scope, matching every other perturbation gate in the "
        "CUDA path).",
        "CUDA_Dihedral_Interaction", io::message::error);
    return 1;
  }
  if (sim.param().gamd.gamd) {
    io::messages.add(
        "CUDA_Dihedral_Interaction does not support GAMD.",
        "CUDA_Dihedral_Interaction", io::message::error);
    return 1;
  }
  if (sim.param().print.monitor_dihedrals) {
    io::messages.add(
        "CUDA_Dihedral_Interaction does not support dihedral-angle-"
        "minimum monitoring (print.monitor_dihedrals).",
        "CUDA_Dihedral_Interaction", io::message::error);
    return 1;
  }

  const std::vector<topology::four_body_term_struct> & dihedrals =
      topo.solute().dihedrals();
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

  const std::vector<interaction::dihedral_type_struct> & types = topo.dihedral_types();
  m_K.resize(types.size());
  m_cospd.resize(types.size());
  m_pd.resize(types.size());
  m_m.resize(types.size());
  for (unsigned t = 0; t < types.size(); ++t) {
    m_K[t] = static_cast<FPL_TYPE>(types[t].K);
    m_cospd[t] = static_cast<FPL_TYPE>(types[t].cospd);
    m_pd[t] = static_cast<FPL_TYPE>(types[t].pd);
    m_m[t] = types[t].m;
  }

  const unsigned num_energy_groups = static_cast<unsigned>(topo.energy_groups().size());
  m_dihedral_energy.resize(num_energy_groups);
  m_virial.resize(9);

  m_touched_atoms = gpu::build_touched_atoms(dihedrals,
      [](const topology::four_body_term_struct & d) {
        return std::array<unsigned, 4>{d.i, d.j, d.k, d.l};
      });

  m_initialized = true;

  if (!quiet)
    os << "CUDA DIHEDRAL INTERACTION\n"
       << "\tterms: " << m_num_dihedrals << "\n"
       << "END\n";
  return 0;
}

int interaction::CUDA_Dihedral_Interaction::calculate_interactions(
    topology::Topology & topo,
    configuration::Configuration & conf,
    simulation::Simulation & sim) {

  m_timer.start(sim);

  if (!m_initialized) {
    m_timer.stop();
    return 1;
  }

  const unsigned num_energy_groups = static_cast<unsigned>(m_dihedral_energy.size());
  cudaMemset(m_dihedral_energy.data(), 0, sizeof(double) * num_energy_groups);
  cudaMemset(m_virial.data(), 0, sizeof(double) * 9);

  const math::CuVArray::View pos =
      sim.cuda().configuration_view(conf, gpu::MIRROR_POS).current().pos;

  static gpu::cuvector<FPL3_TYPE> force;
  const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
  if (force.size() < num_atoms) force.resize(num_atoms);
  cudaMemset(force.data(), 0, sizeof(FPL3_TYPE) * num_atoms);

  gpu::launch_dihedral(
      pos, m_dihedral_i.data(), m_dihedral_j.data(), m_dihedral_k.data(),
      m_dihedral_l.data(), m_dihedral_type.data(),
      m_K.data(), m_cospd.data(), m_pd.data(), m_m.data(),
      m_atom_energy_group.data(), m_num_dihedrals,
      conf.boundary_type, conf.current().box,
      force.data(), m_dihedral_energy.data(), m_virial.data());

  cudaDeviceSynchronize();

  gpu::accumulate_sparse_forces(conf, force.data(), m_touched_atoms);

  for (unsigned g = 0; g < num_energy_groups; ++g) {
    conf.current().energies.dihedral_energy[g] += m_dihedral_energy[g];
  }

  for (unsigned b = 0; b < 3; ++b) {
    for (unsigned a = 0; a < 3; ++a) {
      conf.current().virial_tensor(b, a) += m_virial[b * 3 + a];
    }
  }

  m_timer.stop();
  return 0;
}
