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
 * @file precision_cast_kernels.cu
 * See precision_cast_kernels.h.
 */

#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/memory/types.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/memory/precision.h"

#include "precision_cast_kernels.h"

namespace gpu {

__global__ void cast_fpl3_to_double3_kernel(
    const FPL3_TYPE* __restrict__ src, double3* __restrict__ dst,
    unsigned first, unsigned count) {
  const unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= count) return;
  const unsigned a = first + idx;
  const FPL3_TYPE v = src[a];
  dst[a] = make_double3(static_cast<double>(v.x), static_cast<double>(v.y),
                         static_cast<double>(v.z));
}

__global__ void cast_double3_to_fpl3_kernel(
    const double3* __restrict__ src, FPL3_TYPE* __restrict__ dst,
    unsigned first, unsigned count) {
  const unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= count) return;
  const unsigned a = first + idx;
  const double3 v = src[a];
  dst[a] = FPL3_TYPE{static_cast<FPL_TYPE>(v.x), static_cast<FPL_TYPE>(v.y),
                      static_cast<FPL_TYPE>(v.z)};
}

} // namespace gpu

void gpu::launch_cast_fpl3_to_double3(
    const FPL3_TYPE* src, double3* dst, unsigned first, unsigned count,
    cudaStream_t stream) {
  if (count == 0) return;
  const unsigned threads = 128;
  const unsigned blocks = (count + threads - 1) / threads;
  gpu::cast_fpl3_to_double3_kernel<<<blocks, threads, 0, stream>>>(src, dst, first, count);
}

void gpu::launch_cast_double3_to_fpl3(
    const double3* src, FPL3_TYPE* dst, unsigned first, unsigned count,
    cudaStream_t stream) {
  if (count == 0) return;
  const unsigned threads = 128;
  const unsigned blocks = (count + threads - 1) / threads;
  gpu::cast_double3_to_fpl3_kernel<<<blocks, threads, 0, stream>>>(src, dst, first, count);
}
