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
 * @file constraint_force_publish_kernels.h
 * Publishes a constraint algorithm's own private, per-term constraint_
 * force buffer into the shared GPU mirror's constraint_force field
 * (gpu::MIRROR_CONSTRAINT_FORCE), entirely on-device -- replaces the
 * per-algorithm CPU merge loop (`conf.old().constraint_force(i) +=
 * ...`) every constraint algorithm previously used, which required a
 * host round trip every step. Safe as a plain write (`dst[i] = src[i]
 * * scale`, not atomicAdd): solute vs solvent constraint algorithms
 * write disjoint atom ranges (exactly one of SHAKE/LINCS is active for
 * solute, one of M_SHAKE/SETTLE for solvent), and the mirror's
 * constraint_force is zeroed once per step (CudaManager::
 * zero_mirror_force()) before any of them run.
 *
 * Two shapes, matching the two atom-set shapes constraint algorithms
 * already have (same rationale as velocity_from_delta_kernels.h):
 * `_range` for a contiguous solvent block (M-SHAKE/SETTLE), `_indices`
 * for an explicit, possibly-scattered atom list (SHAKE's
 * constrained_atoms(), covering both solute and solvent atoms it
 * constrains).
 *
 * `scale` exists because each kernel's own private constraint_force
 * convention differs: SETTLE's is already fully scaled (dt2_i baked in
 * at computation time, settle_kernels.cu) so scale=1; SHAKE/M-SHAKE
 * compute a raw, unscaled sum and need scale=1/dt^2 applied here,
 * matching what their own CPU merge loop used to do.
 */

#pragma once

#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/memory/precision.h"

namespace gpu {

  void launch_publish_constraint_force_range(
      FPH3_TYPE* dst, const FPH3_TYPE* src, unsigned first, unsigned count,
      double scale, cudaStream_t stream = 0);

  void launch_publish_constraint_force_indices(
      FPH3_TYPE* dst, const FPH3_TYPE* src, const unsigned* indices,
      unsigned num_indices, double scale, cudaStream_t stream = 0);

  /**
   * Same as launch_publish_constraint_force_range(), but for a source
   * buffer that's unconditionally `double3` regardless of FP_PRECISION
   * (CUDA_Settle's own m_constraint_force -- settle_kernels.cu
   * deliberately stays full double precision, see its own doc
   * comment). Under FP_PRECISION 2/3 (FPH_TYPE == double) this is the
   * same type as the range/indices overloads above; under
   * FP_PRECISION 1 (FPH_TYPE == float) it isn't, hence the separate
   * overload rather than relying on an implicit conversion.
   */
  void launch_publish_constraint_force_range_from_double3(
      FPH3_TYPE* dst, const double3* src, unsigned first, unsigned count,
      double scale, cudaStream_t stream = 0);

} // namespace gpu
