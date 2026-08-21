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
 * @file cuda_shake.cc
 * GPU-native SHAKE constraint algorithm. See cuda_shake.h. Plain C++,
 * no kernel syntax; the real __global__ kernel lives in
 * gpu/cuda/algorithm/constraints/shake_kernels.cu, compiled into
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
#include "cuda_shake.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

int algorithm::CUDA_Shake::init(
    topology::Topology & topo,
    configuration::Configuration & conf,
    simulation::Simulation & sim,
    std::ostream & os,
    bool quiet) {

  if (sim.mpi_enabled()) {
    io::messages.add(
        "CUDA_Shake does not support MPI.",
        "CUDA_Shake", io::message::error);
    return 1;
  }
  if (sim.param().angrest.angrest == simulation::angle_constr ||
      sim.param().dihrest.dihrest == simulation::dihedral_constr) {
    io::messages.add(
        "CUDA_Shake does not support angle/dihedral restraint constraints.",
        "CUDA_Shake", io::message::error);
    return 1;
  }
  if (sim.param().start.shake_pos) {
    io::messages.add(
        "CUDA_Shake does not support shaking initial positions "
        "(start.shake_pos) -- a startup-only cost, use CPU Shake if needed.",
        "CUDA_Shake", io::message::error);
    return 1;
  }

  m_solute_active = topo.solute().distance_constraints().size() &&
      sim.param().constraint.solute.algorithm == simulation::constr_shake &&
      sim.param().constraint.ntc > 1;

  if (!quiet) {
    os << "CUDA_SHAKE\n"
       << "\tsolute\t" << (m_solute_active ? "ON" : "OFF") << "\n";
    if (m_solute_active)
      os << "\t\ttolerance = " << m_solute_tolerance << "\n";
    os << "\tsolvent\t"
       << (sim.param().constraint.solvent.algorithm == simulation::constr_shake ? "ON" : "OFF")
       << "\n";
    if (sim.param().constraint.solvent.algorithm == simulation::constr_shake)
      os << "\t\ttolerance = " << m_solvent_tolerance << "\n";
    os << "END\n";
  }

  const std::vector<interaction::bond_type_struct> & bondtypes = topo.bond_types_harm();

  if (m_solute_active) {
    const std::vector<topology::two_body_term_struct> & dc =
        topo.solute().distance_constraints();
    m_solute_constraints.resize(dc.size());
    for (unsigned c = 0; c < dc.size(); ++c) {
      const double r0 = bondtypes[dc[c].type].r0;
      m_solute_constraints[c] = gpu::ShakeConstraint{dc[c].i, dc[c].j, static_cast<FPL_TYPE>(r0 * r0)};
      constrained_atoms().insert(dc[c].i);
      constrained_atoms().insert(dc[c].j);
    }

    const unsigned num_solute_atoms = static_cast<unsigned>(topo.num_solute_atoms());
    m_solute_inv_mass.resize(num_solute_atoms);
    for (unsigned a = 0; a < num_solute_atoms; ++a) {
      m_solute_inv_mass[a] = static_cast<FPL_TYPE>(topo.inverse_mass()(a));
    }
    m_solute_delta.resize(num_solute_atoms);
    m_changed_flag.resize(1);
  }

  // Only build solvent constraint sets -- and thus only ever run
  // apply()'s solvent SHAKE kernel -- when the solvent algorithm is
  // actually SHAKE. Previously this was unconditional: m_solvent_types
  // got populated from the topology regardless of
  // sim.param().constraint.solvent.algorithm, so apply()'s solvent loop
  // (which iterates m_solvent_types with no algorithm check at all) ran
  // SHAKE on the solvent even when NTCS selected a different algorithm
  // (e.g. SETTLE) -- double-constraining the same water molecules via
  // two independent solvers, corrupting geometry until SHAKE itself
  // failed ("vectors orthogonal"). Matches the CPU Shake class's own
  // gating (shake.cc, checked before every solvent-touching block).
  if (sim.param().constraint.solvent.algorithm == simulation::constr_shake) {
    unsigned first_atom = topo.num_solute_atoms();
    m_solvent_types.resize(topo.num_solvents());
    for (unsigned s = 0; s < topo.num_solvents(); ++s) {
      SolventType & st = m_solvent_types[s];
      st.num_atoms_per_molecule = topo.solvent(s).num_atoms();
      st.num_molecules = topo.num_solvent_molecules(s);
      st.first_atom = first_atom;

      if (st.num_atoms_per_molecule > gpu::MAX_SHAKE_ATOMS_PER_MOLECULE) {
        io::messages.add(
            "CUDA_Shake: solvent molecule exceeds MAX_SHAKE_ATOMS_PER_MOLECULE "
            "(gpu/cuda/algorithm/constraints/shake_kernels.h).",
            "CUDA_Shake", io::message::error);
        return 1;
      }

      const std::vector<topology::two_body_term_struct> & dc =
          topo.solvent(s).distance_constraints();
      st.constraints.resize(dc.size());
      for (unsigned c = 0; c < dc.size(); ++c) {
        const double r0 = bondtypes[dc[c].type].r0;
        st.constraints[c] = gpu::ShakeConstraint{dc[c].i, dc[c].j, static_cast<FPL_TYPE>(r0 * r0)};
      }

      st.inv_mass_local.resize(st.num_atoms_per_molecule);
      for (unsigned a = 0; a < st.num_atoms_per_molecule; ++a) {
        st.inv_mass_local[a] = static_cast<FPL_TYPE>(topo.inverse_mass()(first_atom + a));
      }

      first_atom += st.num_atoms_per_molecule * st.num_molecules;
    }

    for (unsigned int i = topo.num_solute_atoms(); i < topo.num_atoms(); ++i) {
      constrained_atoms().insert(i);
    }
  }

  const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
  m_pos.resize(num_atoms);
  m_old_pos.resize(num_atoms);
  m_constraint_force.resize(num_atoms);
  m_virial.resize(9);
  m_error_flag.resize(1);

  m_initialized = true;
  return 0;
}

int algorithm::CUDA_Shake::apply(
    topology::Topology & topo,
    configuration::Configuration & conf,
    simulation::Simulation & sim) {

  m_timer.start(sim);

  if (!m_initialized) {
    m_timer.stop();
    return 1;
  }

  for (std::set<unsigned int>::const_iterator it = constrained_atoms().begin(),
       to = constrained_atoms().end(); it != to; ++it) {
    conf.old().constraint_force(*it) = 0.0;
  }

  const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
  const unsigned num_solute_atoms = static_cast<unsigned>(topo.num_solute_atoms());

  // FPL_TYPE-narrowing bulk upload (float under FP_PRECISION 1/2) --
  // see vec3_convert.h's doc comment for why this can't be a memcpy.
  gpu::vec3_upload_fpl(m_pos.data(), &conf.current().pos(0), num_atoms);
  gpu::vec3_upload_fpl(m_old_pos.data(), &conf.old().pos(0), num_atoms);
  cudaMemset(m_constraint_force.data(), 0, num_atoms * sizeof(FPH3_TYPE));
  cudaMemset(m_virial.data(), 0, 9 * sizeof(double));
  m_error_flag[0] = 0;

  const double dt = sim.time_step_size();
  const double dt2 = dt * dt;
  const FPL_TYPE dt2_fpl = static_cast<FPL_TYPE>(dt2);

  if (m_solute_active) {
    cudaMemset(m_solute_delta.data(), 0, num_solute_atoms * sizeof(FPL3_TYPE));

    unsigned iterations = 0;
    bool converged = false;
    while (!converged) {
      m_changed_flag[0] = 0;
      gpu::launch_shake_solute_round(
          m_pos.data(), m_old_pos.data(),
          m_solute_constraints.data(), static_cast<unsigned>(m_solute_constraints.size()),
          m_solute_inv_mass.data(), static_cast<FPL_TYPE>(m_solute_tolerance),
          conf.boundary_type, conf.current().box, dt2_fpl,
          m_solute_delta.data(), m_constraint_force.data(), m_virial.data(),
          m_changed_flag.data(), m_error_flag.data());
      gpu::launch_shake_solute_apply(m_pos.data(), m_solute_delta.data(), num_solute_atoms);
      cudaDeviceSynchronize();

      if (m_error_flag[0] != 0) break;
      converged = (m_changed_flag[0] == 0);
      if (++iterations > static_cast<unsigned>(m_max_iterations)) {
        m_error_flag[0] = 2;
        break;
      }
    }

    if (m_error_flag[0] != 0) {
      if (m_error_flag[0] == 1) {
        io::messages.add("SHAKE error. vectors orthogonal",
                          "CUDA_Shake::apply", io::message::error);
        std::cout << "SHAKE failure in solute!" << std::endl;
      } else {
        io::messages.add("SHAKE error. too many iterations",
                          "CUDA_Shake::apply", io::message::critical);
      }
      conf.special().shake_failure_occurred = true;
      m_timer.stop();
      return E_SHAKE_FAILURE_SOLUTE;
    }
  }

  for (const SolventType & st : m_solvent_types) {
    if (st.num_molecules == 0 || st.constraints.size() == 0) continue;
    gpu::launch_shake_solvent(
        m_pos.data(), m_old_pos.data(),
        st.constraints.data(), static_cast<unsigned>(st.constraints.size()),
        st.inv_mass_local.data(), st.num_atoms_per_molecule,
        st.first_atom, st.num_molecules,
        static_cast<FPL_TYPE>(m_solvent_tolerance), static_cast<unsigned>(m_max_iterations),
        conf.boundary_type, conf.current().box, dt2_fpl,
        m_constraint_force.data(), m_virial.data(), m_error_flag.data());
  }
  cudaDeviceSynchronize();

  if (m_error_flag[0] != 0) {
    if (m_error_flag[0] == 1) {
      io::messages.add("SHAKE error. vectors orthogonal",
                        "CUDA_Shake::apply", io::message::error);
      std::cout << "SHAKE failure in solvent!" << std::endl;
    } else {
      io::messages.add("SHAKE error. too many iterations",
                        "CUDA_Shake::apply", io::message::critical);
    }
    conf.special().shake_failure_occurred = true;
    m_timer.stop();
    return E_SHAKE_FAILURE_SOLVENT;
  }

  gpu::vec3_download_fpl(&conf.current().pos(0), m_pos.data(), num_atoms);
  // Accumulation (+=), not a plain copy -- stays a per-atom loop, but
  // only over constrained_atoms() (already exactly the atoms this
  // class's term lists reference), not every atom in the system --
  // same rationale as sparse_force_accumulate.h for the bonded terms.
  // The CPU's own convention (dividing by dt2 once, after the raw
  // lambda*ref_r sum) is preserved here.
  for (unsigned int i : constrained_atoms()) {
    conf.old().constraint_force(i) +=
        math::Vec(static_cast<double>(m_constraint_force[i].x),
                  static_cast<double>(m_constraint_force[i].y),
                  static_cast<double>(m_constraint_force[i].z)) / dt2;
  }

  // Matches the CPU's real `V == math::atomic_virial` gate (shake.h's
  // shake_iteration) -- vacuum boundary never contributes (SPLIT_VIRIAL_
  // BOUNDARY hardcodes math::no_virial for vacuum), and molecular_virial
  // gets nothing added here either (same as CPU): SHAKE's virial
  // contribution is atomic-only, corrected to molecular virial elsewhere
  // if requested (Molecular_Virial_Interaction, generic across
  // accelerators).
  if (conf.boundary_type != math::vacuum &&
      sim.param().pcouple.virial == math::atomic_virial) {
    for (unsigned b = 0; b < 3; ++b) {
      for (unsigned a = 0; a < 3; ++a) {
        conf.old().virial_tensor(b, a) += m_virial[b * 3 + a];
      }
    }
  }

  if (!sim.param().stochastic.sd && !sim.param().minimise.ntem &&
      !sim.param().analyze.analyze) {
    for (std::set<unsigned int>::const_iterator it = constrained_atoms().begin(),
         to = constrained_atoms().end(); it != to; ++it) {
      conf.current().vel(*it) = (conf.current().pos(*it) - conf.old().pos(*it)) / dt;
    }
  }

  m_timer.stop();
  return 0;
}
