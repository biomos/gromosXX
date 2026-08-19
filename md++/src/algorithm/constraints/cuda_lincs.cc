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
 * @file cuda_lincs.cc
 * GPU-native LINCS constraint algorithm. See cuda_lincs.h. Plain C++,
 * no kernel syntax; the real __global__ kernels live in
 * gpu/cuda/algorithm/constraints/lincs_kernels.cu, compiled into
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
#include "lincs.h"
#include "cuda_lincs.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

namespace algorithm {

  /**
   * Uploads one group's static (per-call-time-invariant) data: the
   * term list itself (`constraints`) plus the coupling CSR flattened
   * from `lincs.coupled_constr`/`coef` (topology-owned, built by
   * `setup_lincs()` before this is called). Takes `CUDA_Lincs::Group`'s
   * member cuvectors directly (not the private nested type itself) so
   * it can stay a free function without befriending the class -- the
   * caller fills in the remaining (non-CSR) `Group` fields itself.
   */
  static void fill_group(topology::Topology const & topo,
                          topology::Compound::lincs_struct const & lincs,
                          std::vector<topology::two_body_term_struct> const & constr,
                          unsigned first_atom,
                          gpu::cuvector<gpu::LincsConstraint> & constraints,
                          gpu::cuvector<unsigned> & coupled_offset,
                          gpu::cuvector<unsigned> & coupled_index,
                          gpu::cuvector<FPL_TYPE> & coupled_coef,
                          unsigned & num_constr_per_instance) {
    const unsigned num_constr = static_cast<unsigned>(constr.size());
    num_constr_per_instance = num_constr;

    constraints.resize(num_constr);
    for (unsigned i = 0; i < num_constr; ++i) {
      const double r0 = topo.bond_types_harm()[constr[i].type].r0;
      constraints[i] = gpu::LincsConstraint{
          constr[i].i, constr[i].j,
          static_cast<FPL_TYPE>(r0), static_cast<FPL_TYPE>(lincs.sdiag[i]),
          static_cast<FPL_TYPE>(topo.mass()(constr[i].i + first_atom)),
          static_cast<FPL_TYPE>(topo.mass()(constr[i].j + first_atom))};
    }

    // Flatten lincs.coupled_constr[i]/coef[i] (per-constraint
    // std::vector<unsigned>/<double>) into a single CSR triple.
    coupled_offset.resize(num_constr + 1);
    unsigned total = 0;
    for (unsigned i = 0; i < num_constr; ++i) total += static_cast<unsigned>(lincs.coupled_constr[i].size());
    coupled_index.resize(total);
    coupled_coef.resize(total);

    unsigned pos = 0;
    for (unsigned i = 0; i < num_constr; ++i) {
      coupled_offset[i] = pos;
      for (unsigned n = 0; n < lincs.coupled_constr[i].size(); ++n) {
        coupled_index[pos] = lincs.coupled_constr[i][n];
        coupled_coef[pos] = static_cast<FPL_TYPE>(lincs.coef[i][n]);
        ++pos;
      }
    }
    coupled_offset[num_constr] = pos;
  }

} // namespace algorithm

int algorithm::CUDA_Lincs::init(
    topology::Topology & topo,
    configuration::Configuration & conf,
    simulation::Simulation & sim,
    std::ostream & os,
    bool quiet) {

  if (sim.mpi_enabled()) {
    io::messages.add("CUDA_Lincs does not support MPI.",
                      "CUDA_Lincs", io::message::error);
    return 1;
  }
  if (sim.param().start.shake_pos || sim.param().start.shake_vel) {
    io::messages.add("CUDA_Lincs does not support shaking (lincs-ing) initial "
                      "positions/velocities -- a startup-only cost, use CPU Lincs if needed.",
                      "CUDA_Lincs", io::message::error);
    return 1;
  }

  if (!quiet) {
    os << "CUDA_LINCS\n"
       << "\tsolute\t"
       << (sim.param().constraint.solute.algorithm == simulation::constr_lincs ? "ON" : "OFF")
       << "\n\t\torder = " << sim.param().constraint.solute.lincs_order << "\n"
       << "\tsolvent\t"
       << (sim.param().constraint.solvent.algorithm == simulation::constr_lincs ? "ON" : "OFF")
       << "\n\t\torder = " << sim.param().constraint.solvent.lincs_order << "\n";
  }

  if (sim.param().constraint.solute.algorithm == simulation::constr_lincs) {
    for (auto const & c : topo.solute().distance_constraints()) {
      constrained_atoms().insert(c.i);
      constrained_atoms().insert(c.j);
    }
  }
  if (sim.param().constraint.solvent.algorithm == simulation::constr_lincs) {
    for (unsigned i = topo.num_solute_atoms(); i < topo.num_atoms(); ++i) {
      constrained_atoms().insert(i);
    }
  }

  // Populate topo.solute().lincs()/topo.solvent(i).lincs() (shared,
  // topology-owned setup, same helper the CPU class uses).
  algorithm::setup_lincs(topo, topo.solute().lincs(), topo.solute().distance_constraints());
  unsigned first = static_cast<unsigned>(topo.num_solute_atoms());
  for (unsigned i = 0; i < topo.num_solvents(); ++i) {
    if (topo.num_solvent_molecules(i) != 0) {
      algorithm::setup_lincs(topo, topo.solvent(i).lincs(),
                              topo.solvent(i).distance_constraints(), first);
      first += topo.solvent(i).num_atoms();
    }
  }

  m_solute_active = topo.solute().distance_constraints().size() &&
      sim.param().constraint.solute.algorithm == simulation::constr_lincs &&
      sim.param().constraint.ntc > 1;

  // Sizes every per-round scratch buffer to this group's total work
  // size (num_instances * num_constr_per_instance) -- shared by both
  // the solute and solvent-type group setup below.
  auto resize_scratch = [](Group & g) {
    const unsigned total = g.num_constr_per_instance * g.num_instances;
    g.B.resize(total);
    g.rhs_a.resize(total);
    g.rhs_b.resize(total);
    g.sol.resize(total);
  };

  if (m_solute_active) {
    fill_group(topo, topo.solute().lincs(), topo.solute().distance_constraints(), 0,
               m_solute_group.constraints, m_solute_group.coupled_offset,
               m_solute_group.coupled_index, m_solute_group.coupled_coef,
               m_solute_group.num_constr_per_instance);
    m_solute_group.num_instances = 1;
    m_solute_group.first_atom = 0;
    m_solute_group.atom_stride_per_instance = 0;
    m_solute_group.lincs_order = sim.param().constraint.solute.lincs_order;
    resize_scratch(m_solute_group);
  }

  m_solvent_groups.clear();
  if (sim.param().constraint.solvent.algorithm == simulation::constr_lincs &&
      sim.param().system.nsm) {
    unsigned solvent_first = static_cast<unsigned>(topo.num_solute_atoms());
    for (unsigned s = 0; s < topo.num_solvents(); ++s) {
      const unsigned num_molecules = static_cast<unsigned>(topo.num_solvent_molecules(s));
      if (num_molecules != 0 && topo.solvent(s).distance_constraints().size() > 0) {
        Group g;
        fill_group(topo, topo.solvent(s).lincs(), topo.solvent(s).distance_constraints(),
                   solvent_first, g.constraints, g.coupled_offset, g.coupled_index,
                   g.coupled_coef, g.num_constr_per_instance);
        g.num_instances = num_molecules;
        g.first_atom = solvent_first;
        g.atom_stride_per_instance = static_cast<unsigned>(topo.solvent(s).num_atoms());
        g.lincs_order = sim.param().constraint.solvent.lincs_order;
        resize_scratch(g);
        m_solvent_groups.push_back(std::move(g));
      }
      solvent_first += static_cast<unsigned>(topo.solvent(s).num_atoms()) * num_molecules;
    }
  }

  const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
  m_pos.resize(num_atoms);
  m_old_pos.resize(num_atoms);
  m_rotation_count.resize(1);

  m_initialized = true;
  if (!quiet) os << "END\n";
  return 0;
}

void algorithm::CUDA_Lincs::run_group(
    Group & g, FPL3_TYPE* pos, const FPL3_TYPE* old_pos,
    math::boundary_enum boundary, math::Box box) {

  if (g.num_instances == 0 || g.num_constr_per_instance == 0) return;

  // _solve_lincs's `rec` loop (lincs.cc): each round reads the
  // previous round's solution vector from rhs_a/rhs_b, ping-ponged
  // across calls, and writes rhs_out(i) = sum_n coef*dot*rhs_in(n),
  // sol(i) += rhs_out(i). Shared by the initial solve and the
  // rotational-lengthening correction below -- same structure, just a
  // freshly recomputed starting rhs/sol each time.
  auto run_rounds = [&]() {
    FPL_TYPE * rhs_in = g.rhs_a.data();
    FPL_TYPE * rhs_out = g.rhs_b.data();
    for (int r = 0; r < g.lincs_order; ++r) {
      gpu::launch_lincs_round(g.B.data(), g.coupled_offset.data(), g.coupled_index.data(),
                               g.coupled_coef.data(), g.num_constr_per_instance,
                               g.num_instances, rhs_in, rhs_out, g.sol.data());
      std::swap(rhs_in, rhs_out);
    }
  };

  gpu::launch_lincs_compute_b(old_pos, g.constraints.data(), g.num_constr_per_instance,
                               g.num_instances, g.first_atom, g.atom_stride_per_instance,
                               boundary, box, g.B.data());
  gpu::launch_lincs_init_rhs(pos, g.constraints.data(), g.num_constr_per_instance,
                              g.num_instances, g.first_atom, g.atom_stride_per_instance,
                              boundary, box, g.B.data(), g.rhs_a.data(), g.sol.data());
  run_rounds();
  gpu::launch_lincs_apply(pos, g.constraints.data(), g.B.data(), g.sol.data(),
                           g.num_constr_per_instance, g.num_instances,
                           g.first_atom, g.atom_stride_per_instance);

  // rotational-lengthening correction, second pass -- fresh rhs/sol
  // from the just-updated positions, same round structure again.
  gpu::launch_lincs_rotation_rhs(pos, g.constraints.data(), g.num_constr_per_instance,
                                   g.num_instances, g.first_atom, g.atom_stride_per_instance,
                                   boundary, box, g.rhs_a.data(), g.sol.data(),
                                   m_rotation_count.data());
  run_rounds();
  gpu::launch_lincs_apply(pos, g.constraints.data(), g.B.data(), g.sol.data(),
                           g.num_constr_per_instance, g.num_instances,
                           g.first_atom, g.atom_stride_per_instance);
}

int algorithm::CUDA_Lincs::apply(
    topology::Topology & topo,
    configuration::Configuration & conf,
    simulation::Simulation & sim) {

  m_timer.start(sim);

  if (!m_initialized) {
    m_timer.stop();
    return 1;
  }

  const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
  // FPL_TYPE-narrowing bulk upload (float under FP_PRECISION 1/2) --
  // see vec3_convert.h's doc comment for why this can't be a memcpy.
  gpu::vec3_upload_fpl(m_pos.data(), &conf.current().pos(0), num_atoms);
  gpu::vec3_upload_fpl(m_old_pos.data(), &conf.old().pos(0), num_atoms);
  m_rotation_count[0] = 0;

  if (m_solute_active) {
    run_group(m_solute_group, m_pos.data(), m_old_pos.data(), conf.boundary_type, conf.current().box);
  }
  for (Group & g : m_solvent_groups) {
    run_group(g, m_pos.data(), m_old_pos.data(), conf.boundary_type, conf.current().box);
  }
  cudaDeviceSynchronize();

  gpu::vec3_download_fpl(&conf.current().pos(0), m_pos.data(), num_atoms);

  if (m_rotation_count[0] > 0) {
    std::cout << "LINCS:\ttoo much rotation in " << m_rotation_count[0] << " cases!\n";
  }

  if (!sim.param().stochastic.sd && !sim.param().minimise.ntem &&
      !sim.param().analyze.analyze) {
    for (std::set<unsigned int>::const_iterator it = constrained_atoms().begin(),
         to = constrained_atoms().end(); it != to; ++it) {
      conf.current().vel(*it) = (conf.current().pos(*it) - conf.old().pos(*it)) / sim.time_step_size();
    }
  }

  m_timer.stop();
  return 0;
}
