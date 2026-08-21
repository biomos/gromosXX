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
 * @file berendsen_barostat_kernels.cu
 * GPU-native Berendsen barostat position-scaling kernel.
 */

#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/memory/types.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/utils.h"
#include "math/gmath.h"

#include "berendsen_barostat_kernels.h"

namespace {

__global__ void barostat_scale_positions_kernel(
    math::CuVArray::View pos,
    unsigned num_atoms,
    FPL9_TYPE mu)
{
    const unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_atoms) return;

    const FPL3_TYPE x = pos(i);
    pos(i) = FPL3_TYPE{
        mu(0, 0) * x.x + mu(0, 1) * x.y + mu(0, 2) * x.z,
        mu(1, 0) * x.x + mu(1, 1) * x.y + mu(1, 2) * x.z,
        mu(2, 0) * x.x + mu(2, 1) * x.y + mu(2, 2) * x.z
    };
}

} // namespace

void gpu::launch_barostat_scale_positions(math::CuVArray::View pos,
                                           unsigned num_atoms,
                                           const FPL9_TYPE & mu,
                                           cudaStream_t stream) {
    const unsigned threads = 256;
    const unsigned blocks  = (num_atoms + threads - 1) / threads;
    barostat_scale_positions_kernel<<<blocks, threads, 0, stream>>>(pos, num_atoms, mu);
}
