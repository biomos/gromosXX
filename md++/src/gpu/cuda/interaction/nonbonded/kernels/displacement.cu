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
 * @file displacement.cu
 * Per-atom displacement-since-last-candidate-rebuild kernel. See
 * displacement.h for the host-callable entry point;
 * max_displacement_kernel itself is private to this file (only
 * launch_max_displacement, below, instantiates it).
 */

#include "gpu/cuda/cuheader.h"

#include <algorithm>

#include "gpu/cuda/memory/types.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/math/periodicity.h"
#include "gpu/cuda/utils.h"
#include "math/gmath.h"

#include "displacement.h"

namespace {
  constexpr unsigned NUM_THREADS_PER_BLOCK = 256;
  constexpr unsigned MAX_BLOCKS = 64; // grid-stride loop covers the rest
}

namespace gpu {

// One block-level max per block, via warp-shuffle then a small
// shared-memory reduction across warps -- same two-step shape as
// lj_crf_tile_kernel's energy reduction, just with max instead of sum
// and no atomics needed (each block owns a distinct output slot).
template <math::boundary_enum BOUNDARY>
__global__ void max_displacement_kernel(
    math::CuVArray::View current_pos,
    math::CuVArray::View ref_pos,
    unsigned num_atoms,
    gpu::Periodicity<BOUNDARY> periodicity,
    FPL_TYPE* block_max) {

    const unsigned tid    = threadIdx.x;
    const unsigned idx    = blockIdx.x * blockDim.x + tid;
    const unsigned stride = blockDim.x * gridDim.x;

    FPL_TYPE local_max = FPL_TYPE(0);
    for (unsigned i = idx; i < num_atoms; i += stride) {
        const FPL3_TYPE d = periodicity.nearest_image(current_pos(i), ref_pos(i));
        const FPL_TYPE dist = sqrtf(abs2(d));
        if (dist > local_max) local_max = dist;
    }

    for (unsigned offset = 16; offset > 0; offset >>= 1) {
        const FPL_TYPE other = __shfl_down_sync(0xFFFFFFFFu, local_max, offset);
        if (other > local_max) local_max = other;
    }

    __shared__ FPL_TYPE warp_max[32]; // up to 1024 threads/block == 32 warps
    const unsigned warp_id = tid / 32;
    if ((tid % 32) == 0) warp_max[warp_id] = local_max;
    __syncthreads();

    if (tid == 0) {
        FPL_TYPE result = warp_max[0];
        const unsigned num_warps = (blockDim.x + 31) / 32;
        for (unsigned w = 1; w < num_warps; ++w) {
            if (warp_max[w] > result) result = warp_max[w];
        }
        block_max[blockIdx.x] = result;
    }
}

} // namespace gpu

FPL_TYPE gpu::launch_max_displacement(
    math::CuVArray::View current_pos,
    math::CuVArray::View ref_pos,
    unsigned num_atoms,
    math::boundary_enum boundary,
    math::Box box,
    gpu::cuvector<FPL_TYPE> & partial_max) {

    if (num_atoms == 0) return FPL_TYPE(0);

    const unsigned num_blocks = std::min(
        MAX_BLOCKS, (num_atoms + NUM_THREADS_PER_BLOCK - 1) / NUM_THREADS_PER_BLOCK);
    partial_max.resize(num_blocks);

    switch (boundary) {
        case math::vacuum:
            gpu::max_displacement_kernel<math::vacuum><<<num_blocks, NUM_THREADS_PER_BLOCK>>>(
                current_pos, ref_pos, num_atoms,
                gpu::Periodicity<math::vacuum>(box), partial_max.data());
            break;
        case math::rectangular:
            gpu::max_displacement_kernel<math::rectangular><<<num_blocks, NUM_THREADS_PER_BLOCK>>>(
                current_pos, ref_pos, num_atoms,
                gpu::Periodicity<math::rectangular>(box), partial_max.data());
            break;
        case math::triclinic:
            gpu::max_displacement_kernel<math::triclinic><<<num_blocks, NUM_THREADS_PER_BLOCK>>>(
                current_pos, ref_pos, num_atoms,
                gpu::Periodicity<math::triclinic>(box), partial_max.data());
            break;
        default:
            // Unsupported boundary -- matches CUDA_Pairlist_Algorithm::
            // init()'s vacuum/rectangular-only v1 scope; nothing to
            // launch for other boundaries.
            return FPL_TYPE(0);
    }
    CUDA_CHECK_ERROR("max_displacement_kernel");
    cudaDeviceSynchronize();

    FPL_TYPE result = FPL_TYPE(0);
    for (unsigned b = 0; b < num_blocks; ++b) {
        if (partial_max[b] > result) result = partial_max[b];
    }
    return result;
}
