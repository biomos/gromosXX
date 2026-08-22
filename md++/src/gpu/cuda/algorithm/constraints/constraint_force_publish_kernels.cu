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
 * @file constraint_force_publish_kernels.cu
 * See constraint_force_publish_kernels.h.
 */

#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/memory/types.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/memory/precision.h"

#include "constraint_force_publish_kernels.h"

namespace gpu {

__global__ void publish_constraint_force_range_kernel(
    FPH3_TYPE* __restrict__ dst, const FPH3_TYPE* __restrict__ src,
    unsigned first, unsigned count, double scale) {
  const unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= count) return;
  const unsigned a = first + idx;
  const FPH3_TYPE v = src[a];
  dst[a] = FPH3_TYPE{static_cast<FPH_TYPE>(v.x * scale),
                      static_cast<FPH_TYPE>(v.y * scale),
                      static_cast<FPH_TYPE>(v.z * scale)};
}

__global__ void publish_constraint_force_indices_kernel(
    FPH3_TYPE* __restrict__ dst, const FPH3_TYPE* __restrict__ src,
    const unsigned* __restrict__ indices, unsigned num_indices, double scale) {
  const unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= num_indices) return;
  const unsigned a = indices[idx];
  const FPH3_TYPE v = src[a];
  dst[a] = FPH3_TYPE{static_cast<FPH_TYPE>(v.x * scale),
                      static_cast<FPH_TYPE>(v.y * scale),
                      static_cast<FPH_TYPE>(v.z * scale)};
}

__global__ void publish_constraint_force_range_from_double3_kernel(
    FPH3_TYPE* __restrict__ dst, const double3* __restrict__ src,
    unsigned first, unsigned count, double scale) {
  const unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= count) return;
  const unsigned a = first + idx;
  const double3 v = src[a];
  dst[a] = FPH3_TYPE{static_cast<FPH_TYPE>(v.x * scale),
                      static_cast<FPH_TYPE>(v.y * scale),
                      static_cast<FPH_TYPE>(v.z * scale)};
}

} // namespace gpu

void gpu::launch_publish_constraint_force_range_from_double3(
    FPH3_TYPE* dst, const double3* src, unsigned first, unsigned count,
    double scale, cudaStream_t stream) {
  if (count == 0) return;
  const unsigned threads = 128;
  const unsigned blocks = (count + threads - 1) / threads;
  gpu::publish_constraint_force_range_from_double3_kernel<<<blocks, threads, 0, stream>>>(
      dst, src, first, count, scale);
}

void gpu::launch_publish_constraint_force_range(
    FPH3_TYPE* dst, const FPH3_TYPE* src, unsigned first, unsigned count,
    double scale, cudaStream_t stream) {
  if (count == 0) return;
  const unsigned threads = 128;
  const unsigned blocks = (count + threads - 1) / threads;
  gpu::publish_constraint_force_range_kernel<<<blocks, threads, 0, stream>>>(
      dst, src, first, count, scale);
}

void gpu::launch_publish_constraint_force_indices(
    FPH3_TYPE* dst, const FPH3_TYPE* src, const unsigned* indices,
    unsigned num_indices, double scale, cudaStream_t stream) {
  if (num_indices == 0) return;
  const unsigned threads = 128;
  const unsigned blocks = (num_indices + threads - 1) / threads;
  gpu::publish_constraint_force_indices_kernel<<<blocks, threads, 0, stream>>>(
      dst, src, indices, num_indices, scale);
}
