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
    gpu::TileVecT<gpu::Interaction_Tile>::View candidates) {

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
    bool, gpu::Periodicity<math::vacuum>, FPL_TYPE, gpu::TileVecT<gpu::Interaction_Tile>::View);

template __global__ void gpu::find_block_candidates_kernel<math::rectangular>(
    unsigned, unsigned,
    const FPL3_TYPE*, const FPL_TYPE*,
    const FPL3_TYPE*, const FPL_TYPE*,
    bool, gpu::Periodicity<math::rectangular>, FPL_TYPE, gpu::TileVecT<gpu::Interaction_Tile>::View);

// SPLIT_BOUNDARY (util/template_split.h) routes both truncoct and triclinic
// through the triclinic template instantiation -- needed to link even
// though CUDA_Pairlist_Algorithm::init() hard-errors before this case is
// ever reached at runtime (TILE_PAIRLIST_DESIGN.md §4.1: vacuum +
// rectangular only for v1).
template __global__ void gpu::find_block_candidates_kernel<math::triclinic>(
    unsigned, unsigned,
    const FPL3_TYPE*, const FPL_TYPE*,
    const FPL3_TYPE*, const FPL_TYPE*,
    bool, gpu::Periodicity<math::triclinic>, FPL_TYPE, gpu::TileVecT<gpu::Interaction_Tile>::View);

__global__ void gpu::atom_sort_key_kernel(
    const int* chargegroup_offsets,
    unsigned num_chargegroups,
    unsigned num_atoms,
    const unsigned* cg_sort_key,
    unsigned* atom_sort_key) {

    const unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned stride = blockDim.x * gridDim.x;

    for (unsigned a = idx; a < num_atoms; a += stride) {
        const unsigned cg = gpu::atom_to_chargegroup(chargegroup_offsets, num_chargegroups, a);
        atom_sort_key[a] = cg_sort_key[cg];
    }
}

template <math::boundary_enum BOUNDARY>
__global__ void gpu::classify_tiles_kernel(
    gpu::TileVecT<gpu::Interaction_Tile>::View candidates,
    const unsigned* row_order, unsigned row_count,
    const unsigned* col_other_order, unsigned col_other_count,
    math::CuVArray::View cg_cog,
    gpu::Topology::View topo,
    gpu::Periodicity<BOUNDARY> periodicity,
    FPL_TYPE cutoff_short2, FPL_TYPE cutoff_long2,
    gpu::TileVecT<gpu::Interaction_Tile>::View out_short,
    gpu::TileVecT<gpu::Interaction_Tile>::View out_long) {

    const unsigned tile_idx = blockIdx.x;
    if (tile_idx >= candidates.size()) return;

    const gpu::Interaction_Tile tile = candidates[tile_idx];
    unsigned row_block, col_block;
    bool col_from_b;
    gpu::unpack_block_index(tile.index, row_block, col_block, col_from_b);

    const unsigned r = threadIdx.y;
    const unsigned c = threadIdx.x;

    const unsigned* col_order = col_from_b ? col_other_order : row_order;
    const unsigned  col_count = col_from_b ? col_other_count : row_count;
    // Only a genuine self-pair tile (row/col share the same order array)
    // can contain the same atom twice or a duplicate unordered pair.
    const bool is_diagonal = (!col_from_b) && (row_block == col_block);

    const bool row_valid = (row_block * gpu::BLOCK_SIZE + r) < row_count;
    const bool col_valid = (col_block * gpu::BLOCK_SIZE + c) < col_count;
    const bool skip_dup  = is_diagonal && (r >= c); // avoid double count + self pair

    bool hit_short = false;
    bool hit_long  = false;

    if (row_valid && col_valid && !skip_dup) {
        const unsigned a1 = row_order[row_block * gpu::BLOCK_SIZE + r];
        const unsigned a2 = col_order[col_block * gpu::BLOCK_SIZE + c];

        const unsigned cg1 = gpu::atom_to_chargegroup(topo.chargegroup, topo.num_chargegroups, a1);
        const unsigned cg2 = gpu::atom_to_chargegroup(topo.chargegroup, topo.num_chargegroups, a2);

        if (cg1 != cg2) { // never emit intramolecular (same-chargegroup) pairs
            const unsigned i = a1 < a2 ? a1 : a2;
            const unsigned j = a1 < a2 ? a2 : a1;
            if (!topo.is_excluded(i, j)) {
                const FPL3_TYPE d = periodicity.nearest_image(cg_cog(cg1), cg_cog(cg2));
                const FPL_TYPE dist2 = abs2(d);
                if (dist2 < cutoff_short2) {
                    hit_short = true;
                } else if (dist2 < cutoff_long2) {
                    hit_long = true;
                }
            }
        }
    }

    // blockDim.x == 32 == warp size, threadIdx.y constant per warp -- each
    // warp is exactly one tile row, so __ballot_sync builds that row's
    // mask word with no atomics needed.
    const unsigned short_bits = __ballot_sync(0xFFFFFFFFu, hit_short);
    const unsigned long_bits  = __ballot_sync(0xFFFFFFFFu, hit_long);

    __shared__ unsigned s_mask_short[gpu::BLOCK_SIZE];
    __shared__ unsigned s_mask_long[gpu::BLOCK_SIZE];
    if (c == 0) {
        s_mask_short[r] = short_bits;
        s_mask_long[r]  = long_bits;
    }
    __syncthreads();

    if (r == 0 && c == 0) {
        bool any_short = false, any_long = false;
        for (unsigned k = 0; k < gpu::BLOCK_SIZE; ++k) {
            if (s_mask_short[k]) any_short = true;
            if (s_mask_long[k])  any_long  = true;
        }
        if (any_short) {
            gpu::Interaction_Tile t;
            t.index = tile.index;
            for (unsigned k = 0; k < gpu::BLOCK_SIZE; ++k) t.mask[k] = s_mask_short[k];
            out_short.push_back(t);
        }
        if (any_long) {
            gpu::Interaction_Tile t;
            t.index = tile.index;
            for (unsigned k = 0; k < gpu::BLOCK_SIZE; ++k) t.mask[k] = s_mask_long[k];
            out_long.push_back(t);
        }
    }
}

// explicit instantiations to allow linking
template __global__ void gpu::classify_tiles_kernel<math::vacuum>(
    gpu::TileVecT<gpu::Interaction_Tile>::View,
    const unsigned*, unsigned, const unsigned*, unsigned,
    math::CuVArray::View, gpu::Topology::View,
    gpu::Periodicity<math::vacuum>, FPL_TYPE, FPL_TYPE,
    gpu::TileVecT<gpu::Interaction_Tile>::View, gpu::TileVecT<gpu::Interaction_Tile>::View);

template __global__ void gpu::classify_tiles_kernel<math::rectangular>(
    gpu::TileVecT<gpu::Interaction_Tile>::View,
    const unsigned*, unsigned, const unsigned*, unsigned,
    math::CuVArray::View, gpu::Topology::View,
    gpu::Periodicity<math::rectangular>, FPL_TYPE, FPL_TYPE,
    gpu::TileVecT<gpu::Interaction_Tile>::View, gpu::TileVecT<gpu::Interaction_Tile>::View);

template __global__ void gpu::classify_tiles_kernel<math::triclinic>(
    gpu::TileVecT<gpu::Interaction_Tile>::View,
    const unsigned*, unsigned, const unsigned*, unsigned,
    math::CuVArray::View, gpu::Topology::View,
    gpu::Periodicity<math::triclinic>, FPL_TYPE, FPL_TYPE,
    gpu::TileVecT<gpu::Interaction_Tile>::View, gpu::TileVecT<gpu::Interaction_Tile>::View);
