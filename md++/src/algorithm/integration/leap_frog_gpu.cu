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
 * @file leap_frog_gpu.cu
 * GPU backend specialisations of Leap_Frog_Position and Leap_Frog_Velocity.
 *
 * These kernels operate directly on the GPU-resident CuVArray data held in
 * conf's gpu::Configuration, avoiding CPU↔GPU round-trips every step.
 *
 * After both kernels complete the caller must call
 *   conf.copy_from_gpu()   (full sync, for output steps)
 * or leave data on GPU for the next step.
 */

#include "stdheader.h"

#include "algorithm/algorithm.h"
#include "topology/topology.h"
#include "simulation/simulation.h"
#include "configuration/configuration.h"

#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/configuration_struct.h"
#include "gpu/cuda/utils.h"

#include "leap_frog.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE integration

// ─────────────────────────────────────────────────────────────────────────────
// Device kernels
// ─────────────────────────────────────────────────────────────────────────────

/**
 * v(t+dt/2) = v(t-dt/2) + dt/m * F(t)
 * Operates on GPU arrays directly.
 */
static __global__ void leap_frog_vel_kernel(
    FPL3_TYPE*       cur_vel,
    const FPL3_TYPE* old_vel,
    const FPL3_TYPE* force,
    const float*     inv_mass,
    unsigned         num_atoms,
    FPL_TYPE         dt)
{
    unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
    for (; i < num_atoms; i += gridDim.x * blockDim.x) {
        const FPL_TYPE im = static_cast<FPL_TYPE>(inv_mass[i]);
        cur_vel[i].x = old_vel[i].x + force[i].x * dt * im;
        cur_vel[i].y = old_vel[i].y + force[i].y * dt * im;
        cur_vel[i].z = old_vel[i].z + force[i].z * dt * im;
    }
}

/**
 * x(t+dt) = x(t) + dt * v(t+dt/2)
 * Operates on GPU arrays directly.
 */
static __global__ void leap_frog_pos_kernel(
    FPL3_TYPE*       cur_pos,
    const FPL3_TYPE* old_pos,
    const FPL3_TYPE* cur_vel,
    unsigned         num_atoms,
    FPL_TYPE         dt)
{
    unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
    for (; i < num_atoms; i += gridDim.x * blockDim.x) {
        cur_pos[i].x = old_pos[i].x + cur_vel[i].x * dt;
        cur_pos[i].y = old_pos[i].y + cur_vel[i].y * dt;
        cur_pos[i].z = old_pos[i].z + cur_vel[i].z * dt;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Leap_Frog_Velocity<gpuBackend>
// ─────────────────────────────────────────────────────────────────────────────

template <>
int algorithm::Leap_Frog_Velocity<util::gpuBackend>::apply(
    topology::Topology&           topo,
    configuration::Configuration& conf,
    simulation::Simulation&       sim)
{
    this->m_timer.start(sim);

    // Exchange CPU state (keeps CPU book-keeping correct for output etc.)
    conf.exchange_state();
    conf.current().box = conf.old().box;

    const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
    const FPL_TYPE dt        = static_cast<FPL_TYPE>(sim.time_step_size());

    auto cur = conf.get_gpu_raw_ptrs();
    // After exchange_state, "old" GPU state holds what was "current"
    // We need the old velocity and old force to update current velocity.
    // Sync positions/velocities to GPU (exchange_state only swapped CPU ptrs)
    conf.copy_pos_vel_to_gpu();

    auto old_ptrs = conf.m_gpu->old_raw();

    constexpr unsigned BLOCK = 256;
    const unsigned grid = (num_atoms + BLOCK - 1) / BLOCK;

    leap_frog_vel_kernel<<<grid, BLOCK>>>(
        cur.vel,
        old_ptrs.vel,
        old_ptrs.force,
        topo.get_gpu_view().inverse_mass,
        num_atoms,
        dt);
    CUDA_CHECK_ERROR("leap_frog_vel_kernel");

    // Copy updated velocities back to CPU so other algorithms see them
    CUDA_CHECK(cudaDeviceSynchronize());
    {
        auto& gpu_vel = conf.m_gpu->current.vel;
        for (unsigned i = 0; i < num_atoms; ++i) {
            const auto& v = gpu_vel[i];
            conf.current().vel(i) = math::Vec(v.x, v.y, v.z);
        }
    }

    this->m_timer.stop();
    return 0;
}

// ─────────────────────────────────────────────────────────────────────────────
// Leap_Frog_Position<gpuBackend>
// ─────────────────────────────────────────────────────────────────────────────

template <>
int algorithm::Leap_Frog_Position<util::gpuBackend>::apply(
    topology::Topology&           topo,
    configuration::Configuration& conf,
    simulation::Simulation&       sim)
{
    this->m_timer.start(sim);

    const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
    const FPL_TYPE dt        = static_cast<FPL_TYPE>(sim.time_step_size());

    // Ensure GPU arrays reflect current CPU state
    conf.copy_pos_vel_to_gpu();

    auto cur = conf.get_gpu_raw_ptrs();
    auto old_ptrs = conf.m_gpu->old_raw();

    constexpr unsigned BLOCK = 256;
    const unsigned grid = (num_atoms + BLOCK - 1) / BLOCK;

    leap_frog_pos_kernel<<<grid, BLOCK>>>(
        cur.pos,
        old_ptrs.pos,
        cur.vel,
        num_atoms,
        dt);
    CUDA_CHECK_ERROR("leap_frog_pos_kernel");

    // Copy updated positions back to CPU
    CUDA_CHECK(cudaDeviceSynchronize());
    {
        auto& gpu_pos = conf.m_gpu->current.pos;
        for (unsigned i = 0; i < num_atoms; ++i) {
            const auto& p = gpu_pos[i];
            conf.current().pos(i) = math::Vec(p.x, p.y, p.z);
        }
    }

    this->m_timer.stop();
    return 0;
}

// explicit instantiations
template class algorithm::Leap_Frog_Position<util::gpuBackend>;
template class algorithm::Leap_Frog_Velocity<util::gpuBackend>;
