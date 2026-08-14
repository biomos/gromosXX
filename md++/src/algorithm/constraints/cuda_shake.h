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
 * GPU-native SHAKE constraint algorithm (PLAN.md §10 step 15) -- only
 * ever included under USE_CUDA, see create_constraints.cc, the only
 * call site.
 *
 * v1 scope: SOLVENT distance constraints only. Solvent SHAKE is
 * embarrassingly parallel (every molecule solves its own small,
 * independent Gauss-Seidel iteration -- see shake_kernels.h's doc
 * comment), which is why it's the piece ported here. Solute SHAKE's
 * CPU reference (algorithm::Shake::solute(), shake.h) is a single
 * in-place Gauss-Seidel sweep over *every* solute distance constraint
 * with genuine cross-constraint data dependencies within one iteration
 * -- a fundamentally different (and much harder to parallelize
 * correctly) problem, not attempted here. Hard-errored in init() if
 * actually requested (`ntc > 1` and solute distance constraints exist)
 * rather than silently falling back to a wrong/no-op result -- same
 * convention as every other CUDA scope gate in this codebase. Also
 * hard-errored: MPI (`sim.mpi_enabled()`), angle/dihedral restraint
 * constraints (solute-only, need the same cross-dependency solute
 * solve), and `start.shake_pos` (a startup-only one-time cost, not
 * worth this class's added complexity -- use CPU `Shake` if that
 * combination is needed).
 */

#pragma once

#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/algorithm/constraints/shake_kernels.h"

namespace algorithm {

  class CUDA_Shake : public Algorithm {
  public:
    explicit CUDA_Shake(double const solvent_tolerance = 0.000001,
                         int const max_iterations = 1000,
                         std::string const name = "CUDA_Shake")
      : Algorithm(name),
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

    double m_solvent_tolerance;
    int m_max_iterations;
    std::set<unsigned int> m_constrained_atoms;
    std::vector<SolventType> m_solvent_types;

    gpu::cuvector<double3> m_pos;
    gpu::cuvector<double3> m_old_pos;
    gpu::cuvector<double3> m_constraint_force;
    gpu::cuvector<double> m_virial;
    gpu::cuvector<int> m_error_flag;

    bool m_initialized = false;
  };

} // namespace algorithm
