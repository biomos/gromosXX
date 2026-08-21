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

#include "gpu/cuda/manager/cuda_manager.h"
#include "gpu/constraint_error_slots.h"
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
  m_constraint_force.resize(num_atoms);
  m_virial.resize(9);

  m_constr.resize(3);
  for (unsigned i = 0; i < 3; ++i) m_constr[i] = host_constr[i];
  m_factor_dev.resize(9);
  for (unsigned i = 0; i < 9; ++i) m_factor_dev[i] = static_cast<FPL_TYPE>(m_factor[i]);
  m_constr_length2_dev.resize(3);
  for (unsigned i = 0; i < 3; ++i) m_constr_length2_dev[i] = static_cast<FPL_TYPE>(m_constr_length2[i]);
  m_mass_i_dev.resize(3);
  for (unsigned i = 0; i < 3; ++i) m_mass_i_dev[i] = static_cast<FPL_TYPE>(m_mass_i[i]);

  m_constrained_atoms_dev.resize(constrained_atoms().size());
  {
    unsigned idx = 0;
    for (unsigned int a : constrained_atoms()) m_constrained_atoms_dev[idx++] = a;
  }

  if (m_stream == 0) cudaStreamCreate(&m_stream);

  m_initialized = true;

  if (!quiet) {
    os << "CUDA_M_SHAKE\n"
       << "\tsolvent\ttolerance = " << m_tolerance << "\n"
       << "END\n";
  }
  return 0;
}

algorithm::CUDA_M_Shake::~CUDA_M_Shake() {
  if (m_stream) cudaStreamDestroy(m_stream);
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

  // GPU-resident: reads/writes positions through the shared mirror
  // instead of this algorithm's own private upload/download every
  // call (see cuda_m_shake.h's doc comment for why constraint_force/
  // virial_tensor still use the older private-buffer path for now).
  gpu::Configuration::View view = sim.cuda().configuration_view(conf, gpu::MIRROR_POS, m_stream);

  const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
  const unsigned num_solvent_atoms = num_atoms - m_first_atom;

  cudaMemsetAsync(m_constraint_force.data() + m_first_atom, 0,
                   num_solvent_atoms * sizeof(FPH3_TYPE), m_stream);
  cudaMemsetAsync(m_virial.data(), 0, 9 * sizeof(double), m_stream);

  const double dt = sim.time_step_size();
  const double dt2i = 1.0 / (dt * dt);
  const bool do_virial = sim.param().pcouple.virial == math::atomic_virial;

  int * error_flag = sim.cuda().constraint_error_flag_slot(gpu::ERR_SLOT_M_SHAKE);

  gpu::launch_m_shake_solvent(
      view.current().pos.data(), view.old().pos.data(), m_constr.data(), m_factor_dev.data(),
      m_constr_length2_dev.data(), m_mass_i_dev.data(),
      m_first_atom, m_num_molecules, static_cast<FPL_TYPE>(m_tolerance),
      static_cast<unsigned>(m_max_iterations), static_cast<FPL_TYPE>(dt2i), do_virial,
      m_constraint_force.data(), m_virial.data(), error_flag, m_stream);

  // Vouch for the position we just corrected: no CPU round trip, and
  // gpu_mirror_touches() == 0 keeps Algorithm_Sequence::run()'s
  // default post-apply() invalidation from immediately erasing this.
  sim.cuda().mark_gpu_dirty(conf, gpu::MIRROR_POS, m_stream);

  if (!sim.param().stochastic.sd && !sim.param().minimise.ntem &&
      !sim.param().analyze.analyze) {
    gpu::launch_velocity_from_delta(
        view.current().pos.data(), view.old().pos.data(), view.current().vel.data(),
        m_constrained_atoms_dev.data(), static_cast<unsigned>(m_constrained_atoms_dev.size()),
        static_cast<FPL_TYPE>(1.0 / dt), m_stream);
    sim.cuda().mark_gpu_dirty(conf, gpu::MIRROR_VEL, m_stream);
  }

  // constraint_force/virial_tensor: still the private-buffer path
  // (cuda_m_shake.h's doc comment explains why), so publishing them
  // needs an explicit small sync+download here -- but only on this
  // algorithm's own stream, not a global cudaDeviceSynchronize(), so
  // it doesn't stall whatever else (e.g. CUDA_Lincs) is running
  // concurrently on a different stream.
  cudaStreamSynchronize(m_stream);

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

  m_timer.stop();
  return 0;
}
