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
 * @file remove_com_motion_kernels.cu
 * GPU-native centre-of-mass motion removal kernels.
 *
 * Each reduction kernel has every thread grid-stride-accumulate its own
 * partial sums in registers first, then do exactly one atomicAdd per
 * output element (4, 7, or 12 atomics total per kernel launch -- not
 * per atom). No shared-memory bucketing like the nonbonded tile
 * kernels' energy/virial reductions: this runs at most once per
 * centreofmass.skip_step steps, over a much smaller working set
 * (threads-with-partials, not pairs), so the extra complexity wouldn't
 * pay for itself here.
 */

#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/memory/types.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/utils.h"
#include "math/gmath.h"

#include "remove_com_motion_kernels.h"

namespace {

constexpr unsigned kThreadsPerBlock = 256;

__global__ void com_translation_reduce_kernel(
    math::CuVArray::View vel,
    const float* __restrict__ mass,
    unsigned num_atoms,
    double* sums) {
    double sx = 0, sy = 0, sz = 0, sm = 0;
    for (unsigned i = blockIdx.x * blockDim.x + threadIdx.x; i < num_atoms;
         i += blockDim.x * gridDim.x) {
        const double m = static_cast<double>(mass[i]);
        const FPL3_TYPE v = vel(i);
        sx += m * v.x; sy += m * v.y; sz += m * v.z; sm += m;
    }
    atomicAdd(&sums[0], sx);
    atomicAdd(&sums[1], sy);
    atomicAdd(&sums[2], sz);
    atomicAdd(&sums[3], sm);
}

__global__ void com_translation_apply_kernel(
    math::CuVArray::View vel,
    double com_v_x, double com_v_y, double com_v_z,
    unsigned num_atoms) {
    const unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_atoms) return;
    FPL3_TYPE v = vel(i);
    v.x -= static_cast<FPL_TYPE>(com_v_x);
    v.y -= static_cast<FPL_TYPE>(com_v_y);
    v.z -= static_cast<FPL_TYPE>(com_v_z);
    vel(i) = v;
}

__global__ void com_rotation_reduce_pass1_kernel(
    math::CuVArray::View pos,
    math::CuVArray::View vel,
    const float* __restrict__ mass,
    double dt,
    unsigned num_atoms,
    double* sums) {
    double svx = 0, svy = 0, svz = 0;
    double srx = 0, sry = 0, srz = 0;
    double sm = 0;
    for (unsigned i = blockIdx.x * blockDim.x + threadIdx.x; i < num_atoms;
         i += blockDim.x * gridDim.x) {
        const double m = static_cast<double>(mass[i]);
        const FPL3_TYPE p = pos(i);
        const FPL3_TYPE v = vel(i);
        svx += m * v.x; svy += m * v.y; svz += m * v.z;
        srx += m * p.x - 0.5 * m * v.x * dt;
        sry += m * p.y - 0.5 * m * v.y * dt;
        srz += m * p.z - 0.5 * m * v.z * dt;
        sm  += m;
    }
    atomicAdd(&sums[0], svx);
    atomicAdd(&sums[1], svy);
    atomicAdd(&sums[2], svz);
    atomicAdd(&sums[3], srx);
    atomicAdd(&sums[4], sry);
    atomicAdd(&sums[5], srz);
    atomicAdd(&sums[6], sm);
}

__global__ void com_rotation_reduce_pass2_kernel(
    math::CuVArray::View pos,
    math::CuVArray::View vel,
    const float* __restrict__ mass,
    double dt,
    double com_v_x, double com_v_y, double com_v_z,
    double com_r_x, double com_r_y, double com_r_z,
    unsigned num_atoms,
    double* sums) {
    double Lx = 0, Ly = 0, Lz = 0;
    double I00 = 0, I01 = 0, I02 = 0, I10 = 0, I11 = 0, I12 = 0, I20 = 0, I21 = 0, I22 = 0;
    for (unsigned i = blockIdx.x * blockDim.x + threadIdx.x; i < num_atoms;
         i += blockDim.x * gridDim.x) {
        const double m = static_cast<double>(mass[i]);
        const FPL3_TYPE p = pos(i);
        const FPL3_TYPE v = vel(i);
        const double rx = static_cast<double>(p.x) - 0.5 * dt * static_cast<double>(v.x) - com_r_x;
        const double ry = static_cast<double>(p.y) - 0.5 * dt * static_cast<double>(v.y) - com_r_y;
        const double rz = static_cast<double>(p.z) - 0.5 * dt * static_cast<double>(v.z) - com_r_z;
        const double vx = static_cast<double>(v.x) - com_v_x;
        const double vy = static_cast<double>(v.y) - com_v_y;
        const double vz = static_cast<double>(v.z) - com_v_z;

        // m * cross(r, v - com_v)
        Lx += m * (ry * vz - rz * vy);
        Ly += m * (rz * vx - rx * vz);
        Lz += m * (rx * vy - ry * vx);

        // exact CPU formula (remove_com_motion_cpu.cc): I(0,0)=m*(ry^2+rz^2),
        // I(1,1)=m*(rx^2+rz^2), I(2,2)=m*(rx^2+ry^2), off-diagonals = -m*ri*rj.
        I00 += m * (ry * ry + rz * rz);
        I11 += m * (rx * rx + rz * rz);
        I22 += m * (rx * rx + ry * ry);
        I01 += m * (-rx * ry);
        I10 += m * (-rx * ry);
        I02 += m * (-rx * rz);
        I20 += m * (-rx * rz);
        I12 += m * (-ry * rz);
        I21 += m * (-ry * rz);
    }
    atomicAdd(&sums[0], Lx);
    atomicAdd(&sums[1], Ly);
    atomicAdd(&sums[2], Lz);
    atomicAdd(&sums[3], I00);
    atomicAdd(&sums[4], I01);
    atomicAdd(&sums[5], I02);
    atomicAdd(&sums[6], I10);
    atomicAdd(&sums[7], I11);
    atomicAdd(&sums[8], I12);
    atomicAdd(&sums[9], I20);
    atomicAdd(&sums[10], I21);
    atomicAdd(&sums[11], I22);
}

__global__ void com_rotation_apply_kernel(
    math::CuVArray::View pos,
    math::CuVArray::View vel,
    double dt,
    double com_r_x, double com_r_y, double com_r_z,
    double com_O_x, double com_O_y, double com_O_z,
    unsigned num_atoms) {
    const unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_atoms) return;

    const FPL3_TYPE p = pos(i);
    FPL3_TYPE v = vel(i);
    const double rx = static_cast<double>(p.x) - 0.5 * dt * static_cast<double>(v.x) - com_r_x;
    const double ry = static_cast<double>(p.y) - 0.5 * dt * static_cast<double>(v.y) - com_r_y;
    const double rz = static_cast<double>(p.z) - 0.5 * dt * static_cast<double>(v.z) - com_r_z;

    // cross(com_O, r)
    const double cx = com_O_y * rz - com_O_z * ry;
    const double cy = com_O_z * rx - com_O_x * rz;
    const double cz = com_O_x * ry - com_O_y * rx;

    v.x -= static_cast<FPL_TYPE>(cx);
    v.y -= static_cast<FPL_TYPE>(cy);
    v.z -= static_cast<FPL_TYPE>(cz);
    vel(i) = v;
}

inline unsigned num_blocks_for(unsigned num_atoms) {
    return (num_atoms + kThreadsPerBlock - 1) / kThreadsPerBlock;
}

} // namespace

void gpu::launch_com_translation_reduce(math::CuVArray::View vel,
                                         const float* mass,
                                         unsigned num_atoms,
                                         double* sums) {
    cudaMemset(sums, 0, sizeof(double) * 4);
    const unsigned blocks = num_blocks_for(num_atoms);
    com_translation_reduce_kernel<<<blocks, kThreadsPerBlock>>>(vel, mass, num_atoms, sums);
}

void gpu::launch_com_translation_apply(math::CuVArray::View vel,
                                        double com_v_x, double com_v_y, double com_v_z,
                                        unsigned num_atoms) {
    const unsigned blocks = num_blocks_for(num_atoms);
    com_translation_apply_kernel<<<blocks, kThreadsPerBlock>>>(vel, com_v_x, com_v_y, com_v_z, num_atoms);
}

void gpu::launch_com_rotation_reduce_pass1(math::CuVArray::View pos,
                                           math::CuVArray::View vel,
                                           const float* mass,
                                           double dt,
                                           unsigned num_atoms,
                                           double* sums) {
    cudaMemset(sums, 0, sizeof(double) * 7);
    const unsigned blocks = num_blocks_for(num_atoms);
    com_rotation_reduce_pass1_kernel<<<blocks, kThreadsPerBlock>>>(pos, vel, mass, dt, num_atoms, sums);
}

void gpu::launch_com_rotation_reduce_pass2(math::CuVArray::View pos,
                                           math::CuVArray::View vel,
                                           const float* mass,
                                           double dt,
                                           double com_v_x, double com_v_y, double com_v_z,
                                           double com_r_x, double com_r_y, double com_r_z,
                                           unsigned num_atoms,
                                           double* sums) {
    cudaMemset(sums, 0, sizeof(double) * 12);
    const unsigned blocks = num_blocks_for(num_atoms);
    com_rotation_reduce_pass2_kernel<<<blocks, kThreadsPerBlock>>>(
        pos, vel, mass, dt, com_v_x, com_v_y, com_v_z, com_r_x, com_r_y, com_r_z, num_atoms, sums);
}

void gpu::launch_com_rotation_apply(math::CuVArray::View pos,
                                    math::CuVArray::View vel,
                                    double dt,
                                    double com_r_x, double com_r_y, double com_r_z,
                                    double com_O_x, double com_O_y, double com_O_z,
                                    unsigned num_atoms) {
    const unsigned blocks = num_blocks_for(num_atoms);
    com_rotation_apply_kernel<<<blocks, kThreadsPerBlock>>>(
        pos, vel, dt, com_r_x, com_r_y, com_r_z, com_O_x, com_O_y, com_O_z, num_atoms);
}
