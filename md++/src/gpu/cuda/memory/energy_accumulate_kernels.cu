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
 * @file energy_accumulate_kernels.cu
 * See energy_accumulate_kernels.h.
 */

#include "gpu/cuda/cuheader.h"

#include "energy_accumulate_kernels.h"

namespace gpu {

__global__ void accumulate_energy_kernel(
    double* __restrict__ dst, const double* __restrict__ src, unsigned num_groups) {
  const unsigned i = threadIdx.x;
  if (i >= num_groups) return;
  atomicAdd(&dst[i], src[i]);
}

} // namespace gpu

void gpu::launch_accumulate_energy(double* dst, const double* src, unsigned num_groups, cudaStream_t stream) {
  if (num_groups == 0) return;
  gpu::accumulate_energy_kernel<<<1, num_groups, 0, stream>>>(dst, src, num_groups);
}
