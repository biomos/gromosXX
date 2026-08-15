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
 * @file cuda_shake.h
 * GPU-native SHAKE constraint algorithm (PLAN.md §10 step 15/17) --
 * only ever included under USE_CUDA, see create_constraints.cc, the
 * only call site.
 *
 * Covers both solvent and solute distance constraints, via two
 * genuinely different algorithms (see shake_kernels.h's doc comment
 * for the full derivation of each):
 *  - Solvent: every molecule solves its own small, independent
 *    Gauss-Seidel iteration, bit-comparable to the CPU per molecule
 *    (`launch_shake_solvent`).
 *  - Solute: a Jacobi-style parallel constraint solve (`launch_shake_
 *    solute_round`/`launch_shake_solute_apply`) -- converges to the
 *    same constrained manifold as the CPU's single in-place Gauss-
 *    Seidel sweep (both are standard iterative constraint solvers),
 *    but is **not** bit-comparable to its specific iteration order,
 *    since every constraint in a round reads the same starting
 *    positions instead of the latest updated ones.
 *
 * Hard-errored in init() rather than silently producing a wrong/no-op
 * result: MPI (`sim.mpi_enabled()`), angle/dihedral restraint
 * constraints (solute-only, use a different data structure/solve not
 * ported here), and `start.shake_pos` (a startup-only one-time cost,
 * not worth this class's added complexity -- use CPU `Shake` if that
 * combination is needed).
 */

#pragma once

#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/algorithm/constraints/shake_kernels.h"

namespace algorithm {

  class CUDA_Shake : public Algorithm {
  public:
    explicit CUDA_Shake(double const solute_tolerance = 0.000001,
                         double const solvent_tolerance = 0.000001,
                         int const max_iterations = 1000,
                         std::string const name = "CUDA_Shake")
      : Algorithm(name),
        m_solute_tolerance(solute_tolerance),
        m_solvent_tolerance(solvent_tolerance),
        m_max_iterations(max_iterations) {}

    virtual ~CUDA_Shake() {}

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
     * Per-solvent-type static data, uploaded once in init().
     */
    struct SolventType {
      gpu::cuvector<gpu::ShakeConstraint> constraints;
      gpu::cuvector<double> inv_mass_local;
      unsigned num_atoms_per_molecule = 0;
      unsigned first_atom = 0;
      unsigned num_molecules = 0;
    };

    double m_solute_tolerance;
    double m_solvent_tolerance;
    int m_max_iterations;
    std::set<unsigned int> m_constrained_atoms;
    std::vector<SolventType> m_solvent_types;

    bool m_solute_active = false;
    gpu::cuvector<gpu::ShakeConstraint> m_solute_constraints;
    gpu::cuvector<double> m_solute_inv_mass;
    gpu::cuvector<double3> m_solute_delta;
    gpu::cuvector<int> m_changed_flag;

    gpu::cuvector<double3> m_pos;
    gpu::cuvector<double3> m_old_pos;
    gpu::cuvector<double3> m_constraint_force;
    gpu::cuvector<double> m_virial;
    gpu::cuvector<int> m_error_flag;

    bool m_initialized = false;
  };

} // namespace algorithm
