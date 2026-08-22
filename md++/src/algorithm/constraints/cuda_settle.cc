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
 * @file cuda_settle.cc
 * GPU-native SETTLE constraint algorithm. See cuda_settle.h. Plain
 * C++, no kernel syntax; the real __global__ kernel lives in
 * gpu/cuda/algorithm/constraints/settle_kernels.cu, compiled into
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
#include "gpu/cuda/algorithm/constraints/constraint_force_publish_kernels.h"
#include "gpu/cuda/memory/virial_accumulate_kernels.h"
#include "cuda_settle.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

int algorithm::CUDA_Settle::init(
    topology::Topology & topo,
    configuration::Configuration & conf,
    simulation::Simulation & sim,
    std::ostream & os,
    bool quiet) {

  if (sim.mpi_enabled()) {
    io::messages.add("CUDA_Settle does not support MPI.",
                      "CUDA_Settle", io::message::error);
    return 1;
  }
  if (sim.param().start.shake_pos || sim.param().start.shake_vel) {
    io::messages.add("initial settle-ing is not possible.",
                      "CUDA_Settle", io::message::error);
    return 1;
  }
  if (topo.num_solvents() != 1) {
    io::messages.add("SETTLE does only work if 1 solvent.",
                      "CUDA_Settle", io::message::error);
    return 1;
  }
  if (topo.solvent(0).num_atoms() != 3) {
    io::messages.add("SETTLE does only work with water like molecules (3 atoms).",
                      "CUDA_Settle", io::message::error);
    return 1;
  }
  if (topo.solvent(0).atom(1).mass != topo.solvent(0).atom(2).mass) {
    io::messages.add("SETTLE does only work with water like molecules (wrong masses).",
                      "CUDA_Settle", io::message::error);
    return 1;
  }
  if (topo.solvent(0).distance_constraints().size() != 3) {
    io::messages.add("SETTLE does only work with water like molecules (3 distance constraints).",
                      "CUDA_Settle", io::message::error);
    return 1;
  }
  if (topo.bond_types_harm()[topo.solvent(0).distance_constraint(0).type].r0 !=
      topo.bond_types_harm()[topo.solvent(0).distance_constraint(1).type].r0) {
    io::messages.add("SETTLE does only work with water like molecules (distance constraints wrong).",
                      "CUDA_Settle", io::message::error);
    return 1;
  }

  for (unsigned int i = topo.num_solute_atoms(); i < topo.num_atoms(); ++i) {
    constrained_atoms().insert(i);
  }

  m_mass_O = topo.solvent(0).atom(0).mass;
  m_mass_H = topo.solvent(0).atom(1).mass;
  m_dist_OH = topo.bond_types_harm()[topo.solvent(0).distance_constraint(0).type].r0;
  m_dist_HH = topo.bond_types_harm()[topo.solvent(0).distance_constraint(2).type].r0;
  m_first_atom = static_cast<unsigned>(topo.num_solute_atoms());
  m_num_molecules = static_cast<unsigned>(topo.num_solvent_molecules(0));

  const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
  m_pos.resize(num_atoms);
  m_old_pos.resize(num_atoms);
  m_vel.resize(num_atoms);
  m_constraint_force.resize(num_atoms);
  m_virial.resize(9);
  m_error_flag.resize(1);

  if (m_stream == 0) cudaStreamCreate(&m_stream);

  m_initialized = true;

  if (!quiet) {
    os << "CUDA_SETTLE\n"
       << "\tsolvent\n"
       << "END\n";
  }
  return 0;
}

algorithm::CUDA_Settle::~CUDA_Settle() {
  if (m_stream) cudaStreamDestroy(m_stream);
}

int algorithm::CUDA_Settle::apply(
    topology::Topology & topo,
    configuration::Configuration & conf,
    simulation::Simulation & sim) {

  m_timer.start(sim);

  if (!m_initialized) {
    m_timer.stop();
    return 1;
  }

  if (!(sim.param().system.nsm &&
        sim.param().constraint.solvent.algorithm == simulation::constr_settle)) {
    m_timer.stop();
    return 0;
  }

  const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
  const unsigned num_solvent_atoms = num_atoms - m_first_atom;

  const bool do_velocity = !sim.param().stochastic.sd && !sim.param().minimise.ntem &&
      !sim.param().analyze.analyze;

  // GPU-resident: pos/vel operated on directly through the shared
  // mirror -- no private host upload/download every apply() call (see
  // cuda_settle.h's doc comment). launch_settle() itself stays plain
  // double precision, so the mirror's FPL3_TYPE pos/vel get cast to/
  // from this class's private double3 buffers entirely on-device.
  const unsigned read_fields = do_velocity ? (gpu::MIRROR_POS | gpu::MIRROR_VEL) : gpu::MIRROR_POS;
  gpu::Configuration::View view = sim.cuda().configuration_view(conf, read_fields, m_stream);

  gpu::launch_cast_fpl3_to_double3(view.current().pos.data(), m_pos.data(), m_first_atom, num_solvent_atoms, m_stream);
  gpu::launch_cast_fpl3_to_double3(view.old().pos.data(), m_old_pos.data(), m_first_atom, num_solvent_atoms, m_stream);
  if (do_velocity)
    gpu::launch_cast_fpl3_to_double3(view.current().vel.data(), m_vel.data(), m_first_atom, num_solvent_atoms, m_stream);

  cudaMemsetAsync(m_virial.data(), 0, 9 * sizeof(double), m_stream);
  m_error_flag[0] = 0;

  const double dt_i = 1.0 / sim.time_step_size();

  gpu::launch_settle(
      m_pos.data(), m_old_pos.data(), m_vel.data(), m_first_atom, m_num_molecules,
      m_mass_O, m_mass_H, m_dist_OH, m_dist_HH, dt_i, do_velocity,
      m_constraint_force.data(), m_virial.data(), m_error_flag.data(), m_stream);

  // Cast the corrected positions (and velocities, if computed) back
  // into the mirror -- queued on the same stream right after
  // launch_settle(), so it's covered by the cudaStreamSynchronize()
  // below along with the error-flag check.
  gpu::launch_cast_double3_to_fpl3(m_pos.data(), view.current().pos.data(), m_first_atom, num_solvent_atoms, m_stream);
  if (do_velocity)
    gpu::launch_cast_double3_to_fpl3(m_vel.data(), view.current().vel.data(), m_first_atom, num_solvent_atoms, m_stream);

  cudaStreamSynchronize(m_stream);

  if (m_error_flag[0] != 0) {
    io::messages.add("SETTLE error", "CUDA_Settle", io::message::error);
    std::cout << "SETTLE: exiting with error condition at step " << sim.steps() << std::endl;
    conf.special().shake_failure_occurred = true;
    m_timer.stop();
    return 1;
  }

  // Vouch for what we just corrected: no CPU round trip.
  // gpu_mirror_touches() == 0 keeps Algorithm_Sequence::run()'s default
  // post-apply() invalidation from immediately erasing this.
  sim.cuda().mark_gpu_dirty(conf, gpu::MIRROR_POS, m_stream);
  if (do_velocity) sim.cuda().mark_gpu_dirty(conf, gpu::MIRROR_VEL, m_stream);

  // Publish into the shared mirror's constraint_force (GPU-resident, no
  // CPU round trip) -- a plain write, not atomicAdd: the solvent range
  // is disjoint from whatever solute constraint algorithm is active
  // (SHAKE/LINCS), and the mirror's constraint_force was already
  // zeroed once this step (CudaManager::zero_mirror_force()). scale=1
  // since settle_kernels.cu's constraint_force output is already fully
  // scaled (dt2_i baked in at computation time, unlike SHAKE/M-SHAKE's
  // raw sums).
  gpu::launch_publish_constraint_force_range_from_double3(
      view.old().constraint_force.data(), m_constraint_force.data(),
      m_first_atom, num_solvent_atoms, 1.0, m_stream);
  sim.cuda().mark_gpu_dirty(conf, gpu::MIRROR_CONSTRAINT_FORCE, m_stream);

  // Published via atomicAdd into the shared mirror's virial_tensor
  // (GPU-resident, no CPU round trip) -- writes into old() since
  // exchange_state() already swapped this step, matching
  // constraint_force's own convention above.
  if (sim.param().pcouple.virial == math::atomic_virial) {
    gpu::launch_accumulate_virial9(
        reinterpret_cast<FPH_TYPE*>(view.old().virial_tensor), m_virial.data(), m_stream);
    sim.cuda().mark_gpu_dirty(conf, gpu::MIRROR_VIRIAL, m_stream);
  }

  m_timer.stop();
  return 0;
}
