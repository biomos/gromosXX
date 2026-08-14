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
 * @file lj_crf_tiles.cu
 * LJ + reaction-field force/energy kernel consuming the tile-based GPU
 * pairlist. See lj_crf_tiles.h for the host-callable entry point;
 * `lj_crf_tile_kernel` itself is private to this file (only
 * `launch_lj_crf_tiles`, below, instantiates it), so it needs no public
 * declaration or explicit template instantiation.
 */

#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/memory/types.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/memory/pairlist/tile.h"
#include "gpu/cuda/math/periodicity.h"
#include "gpu/cuda/interaction/nonbonded/pairlist/block_pairlist.h"
#include "gpu/cuda/utils.h"
#include "math/gmath.h"

#include "lj_crf_tiles.h"

namespace gpu {

template <math::boundary_enum BOUNDARY>
__global__ void lj_crf_tile_kernel(
    gpu::TileVecT<gpu::Interaction_Tile>::View tiles,
    const unsigned* row_order, unsigned row_count,
    const unsigned* col_other_order, unsigned col_other_count,
    math::CuVArray::View pos,
    const int* __restrict__ iac,
    const FPL_TYPE* __restrict__ charge,
    const unsigned* __restrict__ atom_energy_group,
    gpu::LJParamView lj,
    gpu::NbSimParams nb,
    gpu::Periodicity<BOUNDARY> periodicity,
    FPL3_TYPE* force,
    double* e_lj_total,
    double* e_crf_total,
    double* virial_total) {

    // Dynamic shared memory, sized (by the launch below) to
    // 2 * num_energy_groups^2 + 9 FPL_TYPEs: [0, G*G) is the LJ bucket
    // matrix, [G*G, 2*G*G) is CRF's, [2*G*G, 2*G*G+9) is the (not
    // energy-group-bucketed) 3x3 atomic virial. G*G == 1 for the common
    // single-energy-group case, degenerating the energy part to exactly
    // the old single-scalar behavior.
    extern __shared__ FPL_TYPE s_energy[];
    const unsigned num_groups  = nb.num_energy_groups;
    const unsigned num_buckets = num_groups * num_groups;
    FPL_TYPE* s_lj     = s_energy;
    FPL_TYPE* s_crf    = s_energy + num_buckets;
    FPL_TYPE* s_virial = s_energy + 2 * num_buckets;

    const unsigned tid = threadIdx.y * blockDim.x + threadIdx.x;
    const unsigned num_stride = blockDim.x * blockDim.y;
    for (unsigned k = tid; k < num_buckets; k += num_stride) {
        s_lj[k]  = 0;
        s_crf[k] = 0;
    }
    for (unsigned k = tid; k < 9; k += num_stride) {
        s_virial[k] = 0;
    }
    __syncthreads();

    const unsigned tile_idx = blockIdx.x;
    if (tile_idx >= tiles.size()) return;

    const gpu::Interaction_Tile tile = tiles[tile_idx];
    unsigned row_block, col_block;
    bool col_from_b;
    gpu::unpack_block_index(tile.index, row_block, col_block, col_from_b);

    const unsigned r = threadIdx.y;
    const unsigned c = threadIdx.x;

    const unsigned* col_order = col_from_b ? col_other_order : row_order;
    const unsigned  col_count = col_from_b ? col_other_count : row_count;

    const bool row_valid = (row_block * gpu::BLOCK_SIZE + r) < row_count;
    const bool col_valid = (col_block * gpu::BLOCK_SIZE + c) < col_count;
    const bool active = row_valid && col_valid && ((tile.mask[r] >> c) & 1u);

    FPL_TYPE e_lj_local  = 0;
    FPL_TYPE e_crf_local = 0;

    if (active) {
        const unsigned a1 = row_order[row_block * gpu::BLOCK_SIZE + r];
        const unsigned a2 = col_order[col_block * gpu::BLOCK_SIZE + c];

        // vector a2 -> a1, matching Nonbonded_Term::lj_crf_interaction's
        // caller convention (nonbonded_innerloop.cc:
        // periodicity.nearest_image(pos(i), pos(j), r)).
        const FPL3_TYPE rvec  = periodicity.nearest_image(pos(a1), pos(a2));
        const FPL_TYPE dist2  = abs2(rvec);
        const FPL_TYPE dist2i = FPL_TYPE(1) / dist2;
        const FPL_TYPE dist6i = dist2i * dist2i * dist2i;

        const FPL2_TYPE ljp = lj.get(iac[a1], iac[a2]);
        const FPL_TYPE c6   = ljp.x;
        const FPL_TYPE c12  = ljp.y;
        const FPL_TYPE c12_dist6i = c12 * dist6i;

        e_lj_local = (c12_dist6i - c6) * dist6i;

        const FPL_TYPE disti = sqrtf(dist2i);
        const FPL_TYPE q_eps = charge[a1] * charge[a2] * nb.four_pi_eps_i;

        e_crf_local = q_eps * (disti - nb.crf_2cut3i * dist2 - nb.crf_cut);

        const FPL_TYPE f = (c12_dist6i + c12_dist6i - c6) * FPL_TYPE(6) * dist6i * dist2i +
                            q_eps * (disti * dist2i + FPL_TYPE(2) * nb.crf_2cut3i);

        const FPL3_TYPE fr = f * rvec;
        atomicAdd(&force[a1].x,  fr.x);
        atomicAdd(&force[a1].y,  fr.y);
        atomicAdd(&force[a1].z,  fr.z);
        atomicAdd(&force[a2].x, -fr.x);
        atomicAdd(&force[a2].y, -fr.y);
        atomicAdd(&force[a2].z, -fr.z);

        // Bucket by energy-group pair, not just tile -- atoms in a tile
        // are sorted by spatial cell, not original index, so there's no
        // shortcut around a per-thread lookup here even though energy
        // groups are contiguous atom-index ranges.
        const unsigned eg_i = atom_energy_group[a1];
        const unsigned eg_j = atom_energy_group[a2];
        const unsigned bucket = eg_i * num_groups + eg_j;
        atomicAdd(&s_lj[bucket],  e_lj_local);
        atomicAdd(&s_crf[bucket], e_crf_local);

        // Atomic virial: virial(b, a) += r(b) * force(a), exact CPU
        // formula (nonbonded_innerloop.cc) -- not energy-group-bucketed,
        // unrolled rather than indexed (FPL3_TYPE has no operator[]).
        atomicAdd(&s_virial[0], rvec.x * fr.x); // (0,0)
        atomicAdd(&s_virial[1], rvec.x * fr.y); // (0,1)
        atomicAdd(&s_virial[2], rvec.x * fr.z); // (0,2)
        atomicAdd(&s_virial[3], rvec.y * fr.x); // (1,0)
        atomicAdd(&s_virial[4], rvec.y * fr.y); // (1,1)
        atomicAdd(&s_virial[5], rvec.y * fr.z); // (1,2)
        atomicAdd(&s_virial[6], rvec.z * fr.x); // (2,0)
        atomicAdd(&s_virial[7], rvec.z * fr.y); // (2,1)
        atomicAdd(&s_virial[8], rvec.z * fr.z); // (2,2)
    }
    __syncthreads();

    // One atomicAdd per bucket (not per pair) into the global totals --
    // degenerates to exactly one atomicAdd per tile, same as before, when
    // num_buckets == 1 (the common single-energy-group case).
    for (unsigned k = tid; k < num_buckets; k += num_stride) {
        atomicAdd(&e_lj_total[k],  static_cast<double>(s_lj[k]));
        atomicAdd(&e_crf_total[k], static_cast<double>(s_crf[k]));
    }
    for (unsigned k = tid; k < 9; k += num_stride) {
        atomicAdd(&virial_total[k], static_cast<double>(s_virial[k]));
    }
}

} // namespace gpu

void gpu::launch_lj_crf_tiles(
    gpu::TileVecT<gpu::Interaction_Tile>::View tiles,
    const unsigned* row_order, unsigned row_count,
    const unsigned* col_other_order, unsigned col_other_count,
    math::CuVArray::View pos,
    const int* iac,
    const FPL_TYPE* charge,
    const unsigned* atom_energy_group,
    gpu::LJParamView lj,
    gpu::NbSimParams nb,
    math::boundary_enum boundary,
    math::Box box,
    FPL3_TYPE* force,
    double* e_lj_total,
    double* e_crf_total,
    double* virial_total,
    cudaStream_t stream) {

    const unsigned num_tiles = tiles.size();
    if (num_tiles == 0) return;

    const dim3 dimBlock(gpu::BLOCK_SIZE, gpu::BLOCK_SIZE);
    const size_t shmem_bytes =
        (2ull * nb.num_energy_groups * nb.num_energy_groups + 9ull) * sizeof(FPL_TYPE);

    switch (boundary) {
        case math::vacuum:
            gpu::lj_crf_tile_kernel<math::vacuum><<<num_tiles, dimBlock, shmem_bytes, stream>>>(
                tiles, row_order, row_count, col_other_order, col_other_count,
                pos, iac, charge, atom_energy_group, lj, nb, gpu::Periodicity<math::vacuum>(box),
                force, e_lj_total, e_crf_total, virial_total);
            break;
        case math::rectangular:
            gpu::lj_crf_tile_kernel<math::rectangular><<<num_tiles, dimBlock, shmem_bytes, stream>>>(
                tiles, row_order, row_count, col_other_order, col_other_count,
                pos, iac, charge, atom_energy_group, lj, nb, gpu::Periodicity<math::rectangular>(box),
                force, e_lj_total, e_crf_total, virial_total);
            break;
        case math::triclinic:
            gpu::lj_crf_tile_kernel<math::triclinic><<<num_tiles, dimBlock, shmem_bytes, stream>>>(
                tiles, row_order, row_count, col_other_order, col_other_count,
                pos, iac, charge, atom_energy_group, lj, nb, gpu::Periodicity<math::triclinic>(box),
                force, e_lj_total, e_crf_total, virial_total);
            break;
        default:
            // Unsupported boundary -- no kernel launch (matches
            // CUDA_Pairlist_Algorithm::init()'s vacuum/rectangular-only
            // v1 scope; there are no tiles for other boundaries anyway).
            break;
    }
    CUDA_CHECK_ERROR("lj_crf_tile_kernel");
}
