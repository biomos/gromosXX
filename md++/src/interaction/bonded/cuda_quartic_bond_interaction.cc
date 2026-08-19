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
 * @file cuda_quartic_bond_interaction.cc
 * GPU-native quartic bond interaction. See cuda_quartic_bond_interaction.h.
 * Plain C++, no kernel syntax; the real __global__ kernel lives in
 * gpu/cuda/interaction/bonded/quartic_bond_kernels.cu, compiled into
 * grocuda.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"
#include "../../interaction/interaction_types.h"

#include "../../gpu/cuda/manager/cuda_manager.h"
#include "gpu/cuda/interaction/bonded/quartic_bond_kernels.h"

#include "cuda_quartic_bond_interaction.h"

#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE bonded

int interaction::CUDA_Quartic_Bond_Interaction::init(
    topology::Topology & topo,
    configuration::Configuration & conf,
    simulation::Simulation & sim,
    std::ostream & os,
    bool quiet) {

  if (sim.param().perturbation.perturbation) {
    io::messages.add(
        "CUDA_Quartic_Bond_Interaction does not support perturbation "
        "(out of scope, matching every other perturbation gate in the "
        "CUDA path).",
        "CUDA_Quartic_Bond_Interaction", io::message::error);
    return 1;
  }
  if (sim.param().gamd.gamd) {
    io::messages.add(
        "CUDA_Quartic_Bond_Interaction does not support GAMD.",
        "CUDA_Quartic_Bond_Interaction", io::message::error);
    return 1;
  }

  const std::vector<topology::two_body_term_struct> & bonds = topo.solute().bonds();
  m_num_bonds = static_cast<unsigned>(bonds.size());

  m_bond_i.resize(m_num_bonds);
  m_bond_j.resize(m_num_bonds);
  m_bond_type.resize(m_num_bonds);
  m_atom_energy_group.resize(m_num_bonds);

  for (unsigned k = 0; k < m_num_bonds; ++k) {
    m_bond_i[k] = static_cast<unsigned>(bonds[k].i);
    m_bond_j[k] = static_cast<unsigned>(bonds[k].j);
    m_bond_type[k] = static_cast<unsigned>(bonds[k].type);
    m_atom_energy_group[k] = topo.atom_energy_group()[bonds[k].i];
  }

  const std::vector<interaction::bond_type_struct> & bondtypes = topo.bond_types_quart();
  m_K.resize(bondtypes.size());
  m_r0.resize(bondtypes.size());
  for (unsigned t = 0; t < bondtypes.size(); ++t) {
    m_K[t] = static_cast<FPL_TYPE>(bondtypes[t].K);
    m_r0[t] = static_cast<FPL_TYPE>(bondtypes[t].r0);
  }

  const unsigned num_energy_groups = static_cast<unsigned>(topo.energy_groups().size());
  m_bond_energy.resize(num_energy_groups);
  m_virial.resize(9);

  if (m_stream == 0) cudaStreamCreate(&m_stream);

  m_initialized = true;

  if (!quiet)
    os << "CUDA QUARTIC BOND INTERACTION\n"
       << "\tterms: " << m_num_bonds << "\n"
       << "END\n";
  return 0;
}

interaction::CUDA_Quartic_Bond_Interaction::~CUDA_Quartic_Bond_Interaction() {
  if (m_stream) cudaStreamDestroy(m_stream);
}

int interaction::CUDA_Quartic_Bond_Interaction::calculate_interactions(
    topology::Topology & topo,
    configuration::Configuration & conf,
    simulation::Simulation & sim) {

  m_timer.start(sim);

  if (!m_initialized) {
    m_timer.stop();
    return 1;
  }

  const unsigned num_energy_groups = static_cast<unsigned>(m_bond_energy.size());
  cudaMemsetAsync(m_bond_energy.data(), 0, sizeof(double) * num_energy_groups, m_stream);
  cudaMemsetAsync(m_virial.data(), 0, sizeof(double) * 9, m_stream);

  // Force written directly into the GPU-resident mirror (zeroed once
  // per step by Forcefield::calculate_interactions()'s sim.cuda().
  // zero_mirror_force(), not here) -- accumulate (+=) via atomicAdd
  // inside the kernel, not overwrite, matching Forcefield's convention.
  gpu::Configuration::View view = sim.cuda().configuration_view(conf, gpu::MIRROR_POS);

  gpu::launch_quartic_bond(
      view.current().pos, m_bond_i.data(), m_bond_j.data(), m_bond_type.data(),
      m_K.data(), m_r0.data(), m_atom_energy_group.data(), m_num_bonds,
      conf.boundary_type, conf.current().box,
      view.current().force.data(), m_bond_energy.data(), m_virial.data(), m_stream);

  sim.cuda().mark_gpu_dirty(conf, gpu::MIRROR_FORCE);

  cudaStreamSynchronize(m_stream);

  for (unsigned g = 0; g < num_energy_groups; ++g) {
    conf.current().energies.bond_energy[g] += m_bond_energy[g];
  }

  for (unsigned b = 0; b < 3; ++b) {
    for (unsigned a = 0; a < 3; ++a) {
      conf.current().virial_tensor(b, a) += m_virial[b * 3 + a];
    }
  }

  m_timer.stop();
  return 0;
}
