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
 * @file cuda_lincs.h
 * GPU-native LINCS constraint algorithm (PLAN.md §10 step 19) -- only
 * ever included under USE_CUDA, see create_constraints.cc, the only
 * call site.
 *
 * Covers both solute (one "group" of `num_instances = 1`, spanning the
 * whole solute system) and solvent (one group per solvent type,
 * `num_instances = num_molecules` of that type, all instances sharing
 * the type's static per-molecule coupling structure) -- see
 * gpu/cuda/algorithm/constraints/lincs_kernels.h for the shared
 * "group" abstraction and why LINCS's own recursion is already
 * Jacobi-shaped (unlike SHAKE), so this port needs no algorithmic
 * reformulation to parallelize, unlike solute SHAKE.
 *
 * v1 scope, hard-errored in init() rather than silently producing a
 * wrong/no-op result: MPI (`sim.mpi_enabled()`) and
 * `start.shake_pos`/`start.shake_vel` (the CPU class's own initial-
 * position/velocity constraining path calls `apply()` recursively
 * during `init()` -- not worth this class's added complexity to
 * replicate for a startup-only cost; use CPU `Lincs` if needed).
 */

#pragma once

#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/algorithm/constraints/lincs_kernels.h"

namespace algorithm {

  class CUDA_Lincs : public Algorithm {
  public:
    CUDA_Lincs() : Algorithm("CUDA_Lincs") {}
    virtual ~CUDA_Lincs() {}

    virtual int init(topology::Topology & topo,
                      configuration::Configuration & conf,
                      simulation::Simulation & sim,
                      std::ostream & os = std::cout,
                      bool quiet = false);

    virtual int apply(topology::Topology & topo,
                       configuration::Configuration & conf,
                       simulation::Simulation & sim);

    std::set<unsigned int> & constrained_atoms() { return m_constrained_atoms; }

  private:
    /**
     * One LINCS "group" (see this file's doc comment): either the
     * whole solute system (`num_instances == 1`) or one solvent type
     * (`num_instances == num_molecules` of that type).
     */
    struct Group {
      gpu::cuvector<gpu::LincsConstraint> constraints;   // per-type, size num_constr_per_instance
      gpu::cuvector<unsigned> coupled_offset;            // per-type CSR, size num_constr_per_instance+1
      gpu::cuvector<unsigned> coupled_index;             // per-type CSR
      gpu::cuvector<FPL_TYPE> coupled_coef;              // per-type CSR
      unsigned num_constr_per_instance = 0;
      unsigned num_instances = 0;
      unsigned first_atom = 0;
      unsigned atom_stride_per_instance = 0;
      int lincs_order = 0;

      // per-call scratch, sized num_instances*num_constr_per_instance
      gpu::cuvector<FPL3_TYPE> B;
      gpu::cuvector<FPL_TYPE> rhs_a;
      gpu::cuvector<FPL_TYPE> rhs_b;
      gpu::cuvector<FPL_TYPE> sol;
    };

    void run_group(Group & g, FPL3_TYPE* pos, const FPL3_TYPE* old_pos,
                   math::boundary_enum boundary, math::Box box);

    std::set<unsigned int> m_constrained_atoms;
    bool m_solute_active = false;
    Group m_solute_group;
    std::vector<Group> m_solvent_groups;

    gpu::cuvector<FPL3_TYPE> m_pos;
    gpu::cuvector<FPL3_TYPE> m_old_pos;
    gpu::cuvector<int> m_rotation_count;

    bool m_initialized = false;
  };

} // namespace algorithm
