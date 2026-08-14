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
 * @file displacement.h
 * Host-callable entry point for the per-atom displacement-since-last-
 * candidate-rebuild kernel (TILE_PAIRLIST_DESIGN.md §10, part 2: the
 * skin-buffer Verlet criterion). Feeds
 * CUDA_Pairlist_Algorithm_Impl::needs_candidate_rebuild(): a candidate
 * rebuild is required once `2 * max_displacement >= skin`, since two
 * atoms could in the worst case have closed the gap between them by
 * that much since the candidates were last built at `cutoff_long + skin`.
 *
 * Same host/device split as lj_crf_tiles.h -- deliberately does not
 * include gpu/cuda/math/periodicity.h or declare the `__global__` kernel
 * itself, since that header isn't host-compilable (bare min/max, only
 * resolve under nvcc). Only this plain launch wrapper is public.
 */

#pragma once

#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "math/gmath.h"
#include "math/box.h"

namespace gpu {

  /**
   * @brief Returns max_i |nearest_image(current_pos[i], ref_pos[i])| --
   * the largest single-atom displacement since `ref_pos` was captured
   * (the last candidate rebuild), measured as a periodicity-wrapped
   * (nearest-image) distance, not a raw coordinate subtraction.
   *
   * Wrapped, not raw, matters: chargegroups get rewrapped into the box
   * every step (prepare_cog()), but `ref_pos` is only refreshed at
   * candidate-rebuild time -- a raw subtraction would see a spurious
   * ~one-box-dimension jump whenever an atom's chargegroup crosses a
   * periodic boundary in between, even though it moved a physically
   * tiny amount. `nearest_image` gives the real physical displacement
   * regardless of which side of the box either snapshot happened to be
   * wrapped to.
   *
   * Synchronous: launches the kernel, synchronizes, and reduces the
   * per-block partial maxima on the host itself (cheap -- this runs at
   * classification cadence, once every `skip_step` steps, not every
   * step, so the extra host round trip is not attempted to be avoided
   * the way the tile kernels' launch sites avoid syncing eagerly).
   *
   * `partial_max` is caller-owned scratch (unified memory), resized by
   * the caller to at least the number of blocks this launches --
   * avoids allocating/freeing it on every call.
   */
  FPL_TYPE launch_max_displacement(math::CuVArray::View current_pos,
                                    math::CuVArray::View ref_pos,
                                    unsigned num_atoms,
                                    math::boundary_enum boundary,
                                    math::Box box,
                                    gpu::cuvector<FPL_TYPE> & partial_max);

}
