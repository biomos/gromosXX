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
 * @file leap_frog_gpu.cc
 * GPU backend specialisations of Leap_Frog_Position and Leap_Frog_Velocity
 * (PLAN.md §10 roadmap step 11). Plain C++ -- groalgorithm has no CUDA
 * language enabled (see remove_com_motion_gpu.cc for the precedent); the
 * real __global__ kernels live in
 * gpu/cuda/algorithm/integration/leap_frog_kernels.cu, compiled into
 * grocuda, and are reached only through the host-safe launch wrappers
 * declared in leap_frog_kernels.h.
 *
 * Velocity leaves its result (new current().vel) resident on the GPU
 * mirror; Position reads that same cached mirror directly (no re-upload,
 * no CPU round trip between the two algorithms) and also computes on the
 * GPU. Only Position's very end syncs back to the CPU-authoritative
 * configuration::Configuration, since no other algorithm in the sequence
 * yet consumes GPU-resident state directly.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../gpu/cuda/manager/cuda_manager.h"
#include "../../gpu/cuda/algorithm/integration/leap_frog_kernels.h"

#include "leap_frog.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE integration

// ─────────────────────────────────────────────────────────────────────────────
// Leap_Frog_Velocity<gpuBackend>
// v(t+dt/2) = v(t-dt/2) + dt/m * F(t)
// ─────────────────────────────────────────────────────────────────────────────

template <>
int algorithm::Leap_Frog_Velocity<util::gpuBackend>::apply(
    topology::Topology&           topo,
    configuration::Configuration& conf,
    simulation::Simulation&       sim)
{
    this->m_timer.start(sim);

    conf.exchange_state();
    conf.current().box = conf.old().box;

    // Full resync: unlike the pairlist path, which only ever needs a
    // cheap pos+vel refresh (force stays GPU-side throughout nonbonded
    // evaluation), the integrator needs old().force fresh too -- and it
    // was last published to the CPU-authoritative conf by the nonbonded
    // interaction's own per-step sync-back. This is the one true
    // resync point per step for the GPU mirror; everything downstream
    // through Leap_Frog_Position stays GPU-resident.
    gpu::Configuration::View view =
        sim.cuda().configuration_view(conf, /*sync_pos_vel=*/false, /*full_resync=*/true);
    const gpu::Topology::View topo_view = sim.cuda().topology_view(topo);

    const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
    const double   dt        = sim.time_step_size();

    gpu::launch_leap_frog_velocity(view.old().vel, view.old().force,
                                    view.current().vel, topo_view.mass,
                                    num_atoms, dt);

    this->m_timer.stop();
    return 0;
}

// ─────────────────────────────────────────────────────────────────────────────
// Leap_Frog_Position<gpuBackend>
// x(t+dt) = x(t) + dt * v(t+dt/2)
// ─────────────────────────────────────────────────────────────────────────────

template <>
int algorithm::Leap_Frog_Position<util::gpuBackend>::apply(
    topology::Topology&           topo,
    configuration::Configuration& conf,
    simulation::Simulation&       sim)
{
    this->m_timer.start(sim);

    if (sim.param().polarise.cos) {
        io::messages.add(
            "Leap_Frog_Position<gpuBackend>: polarisable (COS) charges are "
            "not supported on the GPU integration path.",
            "Leap_Frog_Position", io::message::error);
        this->m_timer.stop();
        return 1;
    }

    const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
    const double   dt        = sim.time_step_size();

    // sync_pos_vel = false, full_resync = false: current().vel was just
    // written by Leap_Frog_Velocity<gpuBackend> on this same cached GPU
    // mirror -- re-syncing from the CPU here would clobber it with the
    // stale value still sitting in conf.current().vel.
    gpu::Configuration::View view =
        sim.cuda().configuration_view(conf, /*sync_pos_vel=*/false, /*full_resync=*/false);

    gpu::launch_leap_frog_position(view.old().pos, view.current().vel,
                                    view.current().pos, num_atoms, dt);

    // The one sync-back per step: publish the GPU-computed pos/vel to
    // the CPU-authoritative conf, since nothing else yet consumes
    // GPU-resident state directly.
    sim.cuda().sync_configuration_from_device(conf);

    this->m_timer.stop();
    return 0;
}

// explicit instantiations for linker
template class algorithm::Leap_Frog_Position<util::gpuBackend>;
template class algorithm::Leap_Frog_Velocity<util::gpuBackend>;
