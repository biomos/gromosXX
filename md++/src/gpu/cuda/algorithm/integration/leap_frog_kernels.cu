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
 * @file leap_frog_kernels.cu
 * GPU-native leap-frog integration kernels (PLAN.md §10 roadmap step 11).
 */

#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/memory/types.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/utils.h"
#include "math/gmath.h"

#include "leap_frog_kernels.h"

namespace {

__global__ void leap_frog_velocity_kernel(
    math::CuVArray::View old_vel,
    math::CuVArray::View old_force,
    math::CuVArray::View new_vel,
    const float* __restrict__ mass,
    unsigned num_atoms,
    FPL_TYPE dt)
{
    const unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_atoms) return;

    const FPL3_TYPE v = old_vel[i];
    const FPL3_TYPE f = old_force[i];
    const FPL_TYPE  m = static_cast<FPL_TYPE>(mass[i]);

    new_vel[i] = FPL3_TYPE{
        v.x + f.x * dt / m,
        v.y + f.y * dt / m,
        v.z + f.z * dt / m
    };
}

__global__ void leap_frog_position_kernel(
    math::CuVArray::View old_pos,
    math::CuVArray::View current_vel,
    math::CuVArray::View new_pos,
    unsigned num_atoms,
    FPL_TYPE dt)
{
    const unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_atoms) return;

    const FPL3_TYPE x = old_pos[i];
    const FPL3_TYPE v = current_vel[i];

    new_pos[i] = FPL3_TYPE{
        x.x + v.x * dt,
        x.y + v.y * dt,
        x.z + v.z * dt
    };
}

} // namespace

void gpu::launch_leap_frog_velocity(math::CuVArray::View old_vel,
                                     math::CuVArray::View old_force,
                                     math::CuVArray::View new_vel,
                                     const float* mass,
                                     unsigned num_atoms,
                                     double dt) {
    const unsigned threads = 256;
    const unsigned blocks  = (num_atoms + threads - 1) / threads;
    leap_frog_velocity_kernel<<<blocks, threads>>>(
        old_vel, old_force, new_vel, mass, num_atoms, static_cast<FPL_TYPE>(dt));
}

void gpu::launch_leap_frog_position(math::CuVArray::View old_pos,
                                     math::CuVArray::View current_vel,
                                     math::CuVArray::View new_pos,
                                     unsigned num_atoms,
                                     double dt) {
    const unsigned threads = 256;
    const unsigned blocks  = (num_atoms + threads - 1) / threads;
    leap_frog_position_kernel<<<blocks, threads>>>(
        old_pos, current_vel, new_pos, num_atoms, static_cast<FPL_TYPE>(dt));
}
