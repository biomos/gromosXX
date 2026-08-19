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
 * @file cuda_m_shake.cc
 * GPU-native M-SHAKE constraint algorithm. See cuda_m_shake.h. Plain
 * C++, no kernel syntax; the real __global__ kernel lives in
 * gpu/cuda/algorithm/constraints/m_shake_kernels.cu, compiled into
 * grocuda.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../interaction/interaction.h"
#include "../../interaction/interaction_types.h"

#include "../../util/error.h"
#include "../../util/debug.h"

#include "gpu/cuda/memory/vec3_convert.h"
#include "cuda_m_shake.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

int algorithm::CUDA_M_Shake::init(
    topology::Topology & topo,
    configuration::Configuration & conf,
    simulation::Simulation & sim,
    std::ostream & os,
    bool quiet) {

  if (sim.mpi_enabled()) {
    io::messages.add("CUDA_M_Shake does not support MPI.",
                      "CUDA_M_Shake", io::message::error);
    return 1;
  }
  if (sim.param().start.shake_pos || sim.param().start.shake_vel) {
    io::messages.add("CUDA_M_Shake does not support shaking initial "
                      "positions/velocities -- use CPU M_Shake if needed.",
                      "CUDA_M_Shake", io::message::error);
    return 1;
  }
  if (topo.num_solvents() != 1) {
    io::messages.add("M-SHAKE is not implemented for multiple solvents.",
                      "CUDA_M_Shake", io::message::error);
    return 1;
  }
  if (topo.solvent(0).distance_constraints().size() != 3) {
    io::messages.add("M-SHAKE is only implemented for 3 solvent constraints.",
                      "CUDA_M_Shake", io::message::error);
    return 1;
  }
  if (topo.solvent(0).atoms().size() != 3) {
    io::messages.add("M-SHAKE is only implemented for 3 solvent atoms.",
                      "CUDA_M_Shake", io::message::error);
    return 1;
  }

  for (unsigned int i = topo.num_solute_atoms(); i < topo.num_atoms(); ++i) {
    constrained_atoms().insert(i);
  }

  m_first_atom = static_cast<unsigned>(topo.num_solute_atoms());
  m_num_molecules = static_cast<unsigned>(topo.num_solvent_molecules(0));

  // Exact port of algorithm::M_Shake::init()'s factor-matrix
  // derivation (m_shake.cc) -- same for every molecule of this
  // (single, required-identical) solvent type, computed once here and
  // reused by every GPU thread every step.
  const std::vector<topology::two_body_term_struct> & dc =
      topo.solvent(0).distance_constraints();

  for (unsigned int i = 0; i < 3; ++i) {
    m_mass_i[i] = 1.0 / topo.mass(topo.num_solute_atoms() + i);
  }

  gpu::MShakeConstraint host_constr[3];
  unsigned k = 0;
  for (auto it_k = dc.begin(), to_k = dc.end(); it_k != to_k; ++it_k, ++k) {
    m_constr_length2[k] = topo.bond_types_harm()[it_k->type].r0 *
                           topo.bond_types_harm()[it_k->type].r0;
    host_constr[k] = gpu::MShakeConstraint{it_k->i, it_k->j};

    unsigned l = 0;
    for (auto it_l = dc.begin(), to_l = dc.end(); it_l != to_l; ++it_l, ++l) {
      int d11 = 0, d12 = 0, d22 = 0, d21 = 0;
      if (it_k->i == it_l->i) {
        d11 = 1;
      } else if (it_k->i == it_l->j) {
        d12 = 1;
      }
      if (it_k->j == it_l->j) {
        d22 = 1;
      } else if (it_k->j == it_l->i) {
        d21 = 1;
      }
      m_factor[3*k+l] = (d11 - d12) * m_mass_i[it_k->i] + (d22 - d21) * m_mass_i[it_k->j];
    }
  }

  const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
  m_pos.resize(num_atoms);
  m_old_pos.resize(num_atoms);
  m_constraint_force.resize(num_atoms);
  m_virial.resize(9);
  m_error_flag.resize(1);

  m_constr.resize(3);
  for (unsigned i = 0; i < 3; ++i) m_constr[i] = host_constr[i];
  m_factor_dev.resize(9);
  for (unsigned i = 0; i < 9; ++i) m_factor_dev[i] = static_cast<FPL_TYPE>(m_factor[i]);
  m_constr_length2_dev.resize(3);
  for (unsigned i = 0; i < 3; ++i) m_constr_length2_dev[i] = static_cast<FPL_TYPE>(m_constr_length2[i]);
  m_mass_i_dev.resize(3);
  for (unsigned i = 0; i < 3; ++i) m_mass_i_dev[i] = static_cast<FPL_TYPE>(m_mass_i[i]);

  m_initialized = true;

  if (!quiet) {
    os << "CUDA_M_SHAKE\n"
       << "\tsolvent\ttolerance = " << m_tolerance << "\n"
       << "END\n";
  }
  return 0;
}

int algorithm::CUDA_M_Shake::apply(
    topology::Topology & topo,
    configuration::Configuration & conf,
    simulation::Simulation & sim) {

  m_timer.start(sim);

  if (!m_initialized) {
    m_timer.stop();
    return 1;
  }

  if (!(sim.param().system.nsm &&
        sim.param().constraint.solvent.algorithm == simulation::constr_m_shake)) {
    m_timer.stop();
    return 0;
  }

  const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
  const unsigned num_solvent_atoms = num_atoms - m_first_atom;

  gpu::vec3_upload_fpl(m_pos.data() + m_first_atom, &conf.current().pos(m_first_atom), num_solvent_atoms);
  gpu::vec3_upload_fpl(m_old_pos.data() + m_first_atom, &conf.old().pos(m_first_atom), num_solvent_atoms);
  cudaMemset(m_constraint_force.data() + m_first_atom, 0, num_solvent_atoms * sizeof(FPL3_TYPE));
  cudaMemset(m_virial.data(), 0, 9 * sizeof(double));
  m_error_flag[0] = 0;

  const double dt = sim.time_step_size();
  const double dt2i = 1.0 / (dt * dt);
  const bool do_virial = sim.param().pcouple.virial == math::atomic_virial;

  gpu::launch_m_shake_solvent(
      m_pos.data(), m_old_pos.data(), m_constr.data(), m_factor_dev.data(),
      m_constr_length2_dev.data(), m_mass_i_dev.data(),
      m_first_atom, m_num_molecules, static_cast<FPL_TYPE>(m_tolerance),
      static_cast<unsigned>(m_max_iterations), static_cast<FPL_TYPE>(dt2i), do_virial,
      m_constraint_force.data(), m_virial.data(), m_error_flag.data());
  cudaDeviceSynchronize();

  if (m_error_flag[0] != 0) {
    if (m_error_flag[0] == 1) {
      io::messages.add("M_SHAKE error. vectors orthogonal",
                        "CUDA_M_Shake::apply", io::message::error);
      std::cout << "M_SHAKE failure in solvent!" << std::endl;
    } else {
      io::messages.add("M_SHAKE error: too many iterations",
                        "CUDA_M_Shake::apply", io::message::critical);
    }
    conf.special().shake_failure_occurred = true;
    m_timer.stop();
    return E_SHAKE_FAILURE_SOLVENT;
  }

  gpu::vec3_download_fpl(&conf.current().pos(m_first_atom), m_pos.data() + m_first_atom, num_solvent_atoms);
  for (unsigned int i : constrained_atoms()) {
    conf.old().constraint_force(i) +=
        math::Vec(static_cast<double>(m_constraint_force[i].x),
                  static_cast<double>(m_constraint_force[i].y),
                  static_cast<double>(m_constraint_force[i].z)) * dt2i;
  }

  if (do_virial) {
    for (unsigned b = 0; b < 3; ++b) {
      for (unsigned a = 0; a < 3; ++a) {
        conf.old().virial_tensor(b, a) += m_virial[b * 3 + a];
      }
    }
  }

  if (!sim.param().stochastic.sd && !sim.param().minimise.ntem &&
      !sim.param().analyze.analyze) {
    const double dti = 1.0 / dt;
    for (unsigned int i : constrained_atoms()) {
      conf.current().vel(i) = (conf.current().pos(i) - conf.old().pos(i)) * dti;
    }
  }

  m_timer.stop();
  return 0;
}
