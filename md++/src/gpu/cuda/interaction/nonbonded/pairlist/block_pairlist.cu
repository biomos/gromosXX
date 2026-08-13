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
 * @file block_pairlist.cu
 * Block bounding-volume + block-pair candidate search kernels.
 */

#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/memory/types.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/memory/pairlist/tile.h"
#include "gpu/cuda/math/periodicity.h"
#include "math/gmath.h"

#include "block_pairlist.h"

__global__ void gpu::compute_block_bounds_kernel(
    const unsigned* order,
    unsigned count,
    unsigned num_blocks,
    math::CuVArray::View cg_cog,
    FPL3_TYPE* block_center,
    FPL_TYPE* block_radius) {

    const unsigned b = blockIdx.x * blockDim.x + threadIdx.x;
    if (b >= num_blocks) return;

    const unsigned begin = b * gpu::BLOCK_SIZE;
    const unsigned end   = min(begin + gpu::BLOCK_SIZE, count);

    FPL3_TYPE center{0, 0, 0};
    for (unsigned k = begin; k < end; ++k) {
        center += cg_cog(order[k]);
    }
    center = center / static_cast<FPL_TYPE>(end - begin);

    FPL_TYPE radius2 = 0;
    for (unsigned k = begin; k < end; ++k) {
        const FPL3_TYPE d = cg_cog(order[k]) - center;
        const FPL_TYPE d2 = abs2(d);
        if (d2 > radius2) radius2 = d2;
    }

    block_center[b] = center;
    block_radius[b] = sqrtf(radius2);
}

template <math::boundary_enum BOUNDARY>
__global__ void gpu::find_block_candidates_kernel(
    unsigned num_blocks_a,
    unsigned num_blocks_b,
    const FPL3_TYPE* block_center_a, const FPL_TYPE* block_radius_a,
    const FPL3_TYPE* block_center_b, const FPL_TYPE* block_radius_b,
    bool self_pairs,
    gpu::Periodicity<BOUNDARY> periodicity,
    FPL_TYPE cutoff,
    gpu::TileVecT<gpu::Interaction_Tile> candidates) {

    const unsigned bi = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned bj = blockIdx.y * blockDim.y + threadIdx.y;

    if (bi >= num_blocks_a || bj >= num_blocks_b) return;
    if (self_pairs && bi > bj) return; // avoid pushing both (i,j) and (j,i)

    const FPL3_TYPE d = periodicity.nearest_image(block_center_a[bi], block_center_b[bj]);
    const FPL_TYPE dist = sqrtf(abs2(d));
    const FPL_TYPE reach = cutoff + block_radius_a[bi] + block_radius_b[bj];

    if (dist > reach) return;

    gpu::Interaction_Tile tile;
    tile.index = gpu::pack_block_index(bi, bj, !self_pairs);
    // mask (exclusions + short/long split) is filled in the classification
    // pass (TILE_PAIRLIST_DESIGN.md §3 step 5), not here -- this kernel's
    // job is exactly "which block-pairs are plausibly close."
    candidates.push_back(tile);
}

// explicit instantiations to allow linking
template __global__ void gpu::find_block_candidates_kernel<math::vacuum>(
    unsigned, unsigned,
    const FPL3_TYPE*, const FPL_TYPE*,
    const FPL3_TYPE*, const FPL_TYPE*,
    bool, gpu::Periodicity<math::vacuum>, FPL_TYPE, gpu::TileVecT<gpu::Interaction_Tile>);

template __global__ void gpu::find_block_candidates_kernel<math::rectangular>(
    unsigned, unsigned,
    const FPL3_TYPE*, const FPL_TYPE*,
    const FPL3_TYPE*, const FPL_TYPE*,
    bool, gpu::Periodicity<math::rectangular>, FPL_TYPE, gpu::TileVecT<gpu::Interaction_Tile>);

// SPLIT_BOUNDARY (util/template_split.h) routes both truncoct and triclinic
// through the triclinic template instantiation -- needed to link even
// though CUDA_Pairlist_Algorithm::init() hard-errors before this case is
// ever reached at runtime (TILE_PAIRLIST_DESIGN.md §4.1: vacuum +
// rectangular only for v1).
template __global__ void gpu::find_block_candidates_kernel<math::triclinic>(
    unsigned, unsigned,
    const FPL3_TYPE*, const FPL_TYPE*,
    const FPL3_TYPE*, const FPL_TYPE*,
    bool, gpu::Periodicity<math::triclinic>, FPL_TYPE, gpu::TileVecT<gpu::Interaction_Tile>);
