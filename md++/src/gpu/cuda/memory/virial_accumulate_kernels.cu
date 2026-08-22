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
 * @file virial_accumulate_kernels.cu
 * See virial_accumulate_kernels.h.
 */

#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/memory/precision.h"

#include "virial_accumulate_kernels.h"

namespace gpu {

__global__ void accumulate_virial9_kernel(
    FPH_TYPE* __restrict__ dst, const double* __restrict__ src) {
  const unsigned i = threadIdx.x;
  if (i >= 9) return;
  atomicAdd(&dst[i], static_cast<FPH_TYPE>(src[i]));
}

} // namespace gpu

void gpu::launch_accumulate_virial9(FPH_TYPE* dst, const double* src, cudaStream_t stream) {
  gpu::accumulate_virial9_kernel<<<1, 9, 0, stream>>>(dst, src);
}
