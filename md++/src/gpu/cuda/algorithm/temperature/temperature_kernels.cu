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
 * @file temperature_kernels.cu
 * GPU-native temperature-group velocity reduction and thermostat
 * velocity-scaling kernels.
 */

#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/memory/types.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/utils.h"
#include "math/gmath.h"

#include "temperature_kernels.h"

namespace {

constexpr unsigned kThreadsPerBlock = 256;

__global__ void group_velocity_reduce_kernel(
    math::CuVArray::View vel,
    const float* __restrict__ mass,
    const unsigned* __restrict__ group_index,
    unsigned num_atoms,
    double* sums) {

    // Direct global atomics into sums[5g+0..4] per temperature group g
    // ([5g+0]=mass, [5g+1..3]=momentum, [5g+4]=self-energy) -- no
    // shared-memory staging. A per-block shared-memory bucket array
    // (the earlier version of this kernel) needed 5*num_groups doubles
    // of *dynamic* shared memory, which is fine for the handful of
    // temperature groups every existing small test topology has, but a
    // real system has one temperature group per (typically rigid,
    // few-atom) solvent molecule -- extended_test/ubiquitin's ~7045
    // groups needs 275KB, blowing well past the ~48KB default dynamic
    // shared memory limit and failing the kernel launch itself
    // (cudaErrorInvalidValue) before a single thread ever runs. Direct
    // atomics have no such ceiling, and contention per group is
    // naturally low precisely because each group is few atoms.
    for (unsigned i = blockIdx.x * blockDim.x + threadIdx.x; i < num_atoms;
         i += blockDim.x * gridDim.x) {
        const double m = static_cast<double>(mass[i]);
        const FPL3_TYPE v = vel(i);
        const unsigned g = group_index[i];
        const unsigned base = 5u * g;
        atomicAdd(&sums[base + 0], m);
        atomicAdd(&sums[base + 1], m * static_cast<double>(v.x));
        atomicAdd(&sums[base + 2], m * static_cast<double>(v.y));
        atomicAdd(&sums[base + 3], m * static_cast<double>(v.z));
        atomicAdd(&sums[base + 4], m * (static_cast<double>(v.x) * v.x +
                                         static_cast<double>(v.y) * v.y +
                                         static_cast<double>(v.z) * v.z));
    }
}

__global__ void thermostat_scale_apply_kernel(
    math::CuVArray::View vel,
    const unsigned* __restrict__ group_index,
    const unsigned* __restrict__ com_bath_of_atom,
    const unsigned* __restrict__ ir_bath_of_atom,
    const FPL3_TYPE* __restrict__ com_v_per_group,
    const double* __restrict__ bath_scale,
    unsigned num_atoms) {
    const unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_atoms) return;

    const FPL3_TYPE com_v = com_v_per_group[group_index[i]];
    const double com_scale = bath_scale[com_bath_of_atom[i]];
    const double ir_scale  = bath_scale[ir_bath_of_atom[i]];

    FPL3_TYPE v = vel(i);
    const double irx = static_cast<double>(v.x) - com_v.x;
    const double iry = static_cast<double>(v.y) - com_v.y;
    const double irz = static_cast<double>(v.z) - com_v.z;

    v.x = static_cast<FPL_TYPE>(com_scale * com_v.x + ir_scale * irx);
    v.y = static_cast<FPL_TYPE>(com_scale * com_v.y + ir_scale * iry);
    v.z = static_cast<FPL_TYPE>(com_scale * com_v.z + ir_scale * irz);
    vel(i) = v;
}

inline unsigned num_blocks_for(unsigned num_atoms) {
    return (num_atoms + kThreadsPerBlock - 1) / kThreadsPerBlock;
}

} // namespace

void gpu::launch_group_velocity_reduce(math::CuVArray::View vel,
                                        const float* mass,
                                        const unsigned* group_index,
                                        unsigned num_atoms,
                                        unsigned num_groups,
                                        double* sums) {
    const size_t sums_bytes = 5ull * num_groups * sizeof(double);
    cudaMemset(sums, 0, sums_bytes);
    const unsigned blocks = num_blocks_for(num_atoms);
    group_velocity_reduce_kernel<<<blocks, kThreadsPerBlock>>>(
        vel, mass, group_index, num_atoms, sums);
}

void gpu::launch_thermostat_scale_apply(math::CuVArray::View vel,
                                         const unsigned* group_index,
                                         const unsigned* com_bath_of_atom,
                                         const unsigned* ir_bath_of_atom,
                                         const FPL3_TYPE* com_v_per_group,
                                         const double* bath_scale,
                                         unsigned num_atoms) {
    const unsigned blocks = num_blocks_for(num_atoms);
    thermostat_scale_apply_kernel<<<blocks, kThreadsPerBlock>>>(
        vel, group_index, com_bath_of_atom, ir_bath_of_atom,
        com_v_per_group, bath_scale, num_atoms);
}
