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

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"
#include "../../interaction/interaction_types.h"

#include "../../gpu/cuda/manager/cuda_manager.h"
#include "gpu/cuda/memory/virial_accumulate_kernels.h"
#include "gpu/cuda/memory/energy_accumulate_kernels.h"
#include "gpu/cuda/interaction/bonded/dihedral_kernels.h"

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

  if (m_stream == 0) cudaStreamCreate(&m_stream);

  m_initialized = true;

  if (!quiet)
    os << "CUDA DIHEDRAL INTERACTION\n"
       << "\tterms: " << m_num_dihedrals << "\n"
       << "END\n";
  return 0;
}

interaction::CUDA_Dihedral_Interaction::~CUDA_Dihedral_Interaction() {
  if (m_stream) cudaStreamDestroy(m_stream);
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
  cudaMemsetAsync(m_dihedral_energy.data(), 0, sizeof(double) * num_energy_groups, m_stream);
  cudaMemsetAsync(m_virial.data(), 0, sizeof(double) * 9, m_stream);

  // Force written directly into the GPU-resident mirror (zeroed once
  // per step by Forcefield::calculate_interactions()'s sim.cuda().
  // zero_mirror_force(), not here) -- no private scratch buffer, no
  // sync, no host readback. Only requesting MIRROR_POS as a read field
  // (never MIRROR_FORCE) is what keeps this from triggering a coarse
  // resync that would clobber force NonBonded/other bonded terms may
  // have already accumulated this step.
  gpu::Configuration::View view = sim.cuda().configuration_view(conf, gpu::MIRROR_POS, m_stream);

  if (!m_energy_registered) {
    sim.cuda().ensure_energy_groups(conf, num_energy_groups);
    m_energy_registered = true;
  }

  gpu::launch_dihedral(
      view.current().pos, m_dihedral_i.data(), m_dihedral_j.data(), m_dihedral_k.data(),
      m_dihedral_l.data(), m_dihedral_type.data(),
      m_K.data(), m_cospd.data(), m_pd.data(), m_m.data(),
      m_atom_energy_group.data(), m_num_dihedrals,
      conf.boundary_type, conf.current().box,
      view.current().force.data(), m_dihedral_energy.data(), m_virial.data(), m_stream);

  sim.cuda().mark_gpu_dirty(conf, gpu::MIRROR_FORCE, m_stream);

  // Publish virial into the shared mirror's virial_tensor via atomicAdd
  // (GPU-resident, no CPU round trip) -- see cuda_angle_interaction.cc
  // for the full rationale (every contributor shares one accumulator,
  // zeroed once per step by CudaManager::zero_mirror_force()).
  gpu::launch_accumulate_virial9(
      reinterpret_cast<FPH_TYPE*>(view.current().virial_tensor), m_virial.data(), m_stream);
  sim.cuda().mark_gpu_dirty(conf, gpu::MIRROR_VIRIAL, m_stream);

  // Same on-device merge for energy -- see energy_accumulate_kernels.h
  // and cuda_angle_interaction.cc's identical comment. No host sync at
  // all in this class anymore.
  const gpu::EnergyMirrorPtrs eptrs = sim.cuda().energy_mirror_ptrs(conf);
  gpu::launch_accumulate_energy(eptrs.dihedral, m_dihedral_energy.data(), num_energy_groups, m_stream);
  sim.cuda().mark_gpu_dirty(conf, gpu::MIRROR_ENERGY, m_stream);

  m_timer.stop();
  return 0;
}
