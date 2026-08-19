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
 * @file velocity_from_delta_kernels.cu
 * See velocity_from_delta_kernels.h.
 */

#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/memory/types.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/memory/precision.h"

#include "velocity_from_delta_kernels.h"

namespace gpu {

__global__ void velocity_from_delta_kernel(
    const FPL3_TYPE* __restrict__ pos,
    const FPL3_TYPE* __restrict__ old_pos,
    FPL3_TYPE* __restrict__ vel,
    const unsigned* __restrict__ atom_indices,
    unsigned num_indices,
    FPL_TYPE dt_i) {

  const unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= num_indices) return;

  const unsigned a = atom_indices[idx];
  vel[a] = (pos[a] - old_pos[a]) * dt_i;
}

} // namespace gpu

void gpu::launch_velocity_from_delta(
    const FPL3_TYPE* pos,
    const FPL3_TYPE* old_pos,
    FPL3_TYPE* vel,
    const unsigned* atom_indices,
    unsigned num_indices,
    FPL_TYPE dt_i,
    cudaStream_t stream) {

  if (num_indices == 0) return;
  const unsigned threads = 128;
  const unsigned blocks = (num_indices + threads - 1) / threads;
  gpu::velocity_from_delta_kernel<<<blocks, threads, 0, stream>>>(
      pos, old_pos, vel, atom_indices, num_indices, dt_i);
}
