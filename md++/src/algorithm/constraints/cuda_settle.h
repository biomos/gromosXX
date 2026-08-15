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
 * @file cuda_settle.h
 * GPU-native SETTLE constraint algorithm (PLAN.md §10 step 18) --
 * only ever included under USE_CUDA, see create_constraints.cc, the
 * only call site.
 *
 * SETTLE is analytical (no iteration), and every molecule's solve is
 * fully independent of every other's -- one GPU thread per molecule
 * doing exactly the CPU's closed-form vector algebra (see
 * gpu/cuda/algorithm/constraints/settle_kernels.h) is both correct and
 * embarrassingly parallel, and *is* bit-comparable to the CPU result
 * per molecule (unlike solute SHAKE's Jacobi scheme).
 *
 * v1 scope, hard-errored in init() rather than silently producing a
 * wrong/no-op result -- same conditions algorithm::Settle::init()
 * itself checks: exactly one solvent type, exactly 3 atoms per
 * molecule, H1/H2 same mass, exactly 3 distance constraints with the
 * two O-H constraints sharing one length. Also hard-errored: MPI
 * (`sim.mpi_enabled()`) and `start.shake_pos`/`start.shake_vel` (the
 * CPU class itself refuses these -- "initial settle-ing is not
 * possible").
 */

#pragma once

#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/algorithm/constraints/settle_kernels.h"

namespace algorithm {

  class CUDA_Settle : public Algorithm {
  public:
    CUDA_Settle() : Algorithm("CUDA_Settle") {}
    virtual ~CUDA_Settle() {}

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
    std::set<unsigned int> m_constrained_atoms;

    double m_mass_O = 0.0;
    double m_mass_H = 0.0;
    double m_dist_OH = 0.0;
    double m_dist_HH = 0.0;
    unsigned m_first_atom = 0;
    unsigned m_num_molecules = 0;

    gpu::cuvector<double3> m_pos;
    gpu::cuvector<double3> m_old_pos;
    gpu::cuvector<double3> m_vel;
    gpu::cuvector<double3> m_constraint_force;
    gpu::cuvector<double> m_virial;
    gpu::cuvector<int> m_error_flag;

    bool m_initialized = false;
  };

} // namespace algorithm
