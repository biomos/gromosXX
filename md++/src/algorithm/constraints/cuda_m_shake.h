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
 * @file cuda_m_shake.h
 * GPU-native M-SHAKE constraint algorithm (solvent only -- solute
 * stays on CUDA_Shake/CUDA_Settle/CUDA_Lincs, matching M_Shake's own
 * CPU-side scope, see create_constraints.cc) -- only ever included
 * under USE_CUDA, see create_constraints.cc, the only call site.
 *
 * M-SHAKE solves all 3 constraints of a molecule simultaneously each
 * iteration via a direct 3x3 matrix inversion (algorithm::M_Shake,
 * m_shake.cc), rather than SHAKE's one-constraint-at-a-time
 * Gauss-Seidel sweep -- converges in far fewer iterations, which
 * matters specifically because profiling CUDA_Shake's solvent path
 * found its cost dominated by per-molecule iteration count under a
 * grid too small to hide each iteration's global-memory latency
 * (KNOWN_ISSUES.md) -- fewer iterations means a shorter unhideable
 * per-thread dependency chain, not just less arithmetic.
 *
 * GPU-resident: reads/writes positions through the shared
 * CudaManager mirror (sim.cuda().configuration_view()/mark_gpu_dirty())
 * instead of its own private upload/download every apply() call, and
 * runs on its own stream with no cudaDeviceSynchronize() at all --
 * status is written into the shared deferred constraint-error-flags
 * buffer (gpu/constraint_error_slots.h) and checked once, centrally,
 * at the end of Algorithm_Sequence::run(). This lets it run genuinely
 * concurrently with CUDA_Lincs's solute pass (disjoint atom ranges,
 * separate streams, no data hazard) instead of serializing through a
 * host round-trip between them.
 *
 * constraint_force publishes into the shared mirror's MIRROR_
 * CONSTRAINT_FORCE field on-device (constraint_force_publish_
 * kernels.h) -- a plain write, since solute/solvent constraint
 * algorithms always write disjoint atom ranges. virial_tensor
 * publishes via atomicAdd (virial_accumulate_kernels.h) instead, since
 * it's a single global accumulator every bonded/nonbonded term also
 * contributes to in the same step -- both are fully GPU-resident now,
 * no CPU round trip in either case.
 *
 * v1 scope, hard-errored in init() rather than silently producing a
 * wrong/no-op result -- same conditions algorithm::M_Shake::init()
 * itself checks: exactly one solvent type, exactly 3 atoms per
 * molecule, exactly 3 distance constraints. Also hard-errored: MPI
 * (`sim.mpi_enabled()`) and `start.shake_pos`/`start.shake_vel` (the
 * CPU class's own initial-shake support isn't ported).
 */

#pragma once

#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/algorithm/constraints/m_shake_kernels.h"
#include "gpu/cuda/algorithm/constraints/velocity_from_delta_kernels.h"

namespace algorithm {

  class CUDA_M_Shake : public Algorithm {
  public:
    CUDA_M_Shake(double const tolerance = 0.000001, int const max_iterations = 1000)
      : Algorithm("CUDA_M_Shake"), m_tolerance(tolerance), m_max_iterations(max_iterations) {}
    virtual ~CUDA_M_Shake();

    virtual int init(topology::Topology & topo,
                      configuration::Configuration & conf,
                      simulation::Simulation & sim,
                      std::ostream & os = std::cout,
                      bool quiet = false);

    virtual int apply(topology::Topology & topo,
                       configuration::Configuration & conf,
                       simulation::Simulation & sim);

    std::set<unsigned int> & constrained_atoms() { return m_constrained_atoms; }

    // Owns its own GPU-mirror freshness: apply() marks MIRROR_POS
    // fresh+dirty itself and must not have that immediately erased by
    // Algorithm_Sequence::run()'s default post-apply() invalidation
    // (see leap_frog_gpu.cc for the same pattern).
    virtual unsigned gpu_mirror_touches() const { return 0; }

  private:
    std::set<unsigned int> m_constrained_atoms;

    double m_tolerance;
    int m_max_iterations;

    // Same for every molecule (single, identical solvent type) --
    // computed once in init(), matches algorithm::M_Shake's own
    // factor/constr_length2/mass_i member variables exactly.
    double m_factor[9] = {0,0,0,0,0,0,0,0,0};      // row-major 3x3
    double m_constr_length2[3] = {0,0,0};
    double m_mass_i[3] = {0,0,0};

    unsigned m_first_atom = 0;
    unsigned m_num_molecules = 0;

    gpu::cuvector<gpu::MShakeConstraint> m_constr;
    gpu::cuvector<FPL_TYPE> m_factor_dev;
    gpu::cuvector<FPL_TYPE> m_constr_length2_dev;
    gpu::cuvector<FPL_TYPE> m_mass_i_dev;
    gpu::cuvector<FPH3_TYPE> m_constraint_force;
    gpu::cuvector<double> m_virial;
    // constrained_atoms(), uploaded once in init(), for the on-device
    // velocity_from_delta kernel -- contiguous here (the whole solvent
    // range) but stored explicitly anyway to share the one kernel with
    // CUDA_Lincs's possibly-non-contiguous solute chain.
    gpu::cuvector<unsigned> m_constrained_atoms_dev;

    // Own stream: lets this run concurrently with CUDA_Lincs's solute
    // pass (or anything else) instead of implicitly serializing on the
    // default stream.
    cudaStream_t m_stream = 0;

    bool m_initialized = false;
  };

} // namespace algorithm
