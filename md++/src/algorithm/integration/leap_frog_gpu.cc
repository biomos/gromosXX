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
    // Mirror's own current/old halves swap at the exact same point the
    // CPU-authoritative conf's do -- without this, view.old().force/
    // .old().vel below would not actually refer to the buffer Forcefield/
    // the previous step's velocity write just filled (see CudaManager::
    // exchange_mirror_state()'s doc comment).
    sim.cuda().exchange_mirror_state(conf);
    conf.current().box = conf.old().box;

    // Own stream: MIRROR_FORCE is typically produced by five different
    // algorithms (four bonded terms + NonBonded, each on their own
    // stream) -- passing this stream to configuration_view() below lets
    // it insert cudaStreamWaitEvent() against all of them instead of
    // relying on legacy-default-stream implicit ordering.
    static cudaStream_t stream = 0;
    if (stream == 0) cudaStreamCreate(&stream);

    // MIRROR_BOX deliberately not requested: launch_leap_frog_velocity()
    // takes no box argument at all (confirmed against leap_frog_kernels.h)
    // -- it was dead weight in this mask, and forced configuration_view()'s
    // coarse full-upload path every step (no per-field upload routine
    // covers BOX), silently masking the fact that view.old().force wasn't
    // actually being kept correct any other way (fixed above).
    gpu::Configuration::View view = sim.cuda().configuration_view(
        conf, gpu::MIRROR_POS | gpu::MIRROR_VEL | gpu::MIRROR_FORCE, stream);
    const gpu::Topology::View topo_view = sim.cuda().topology_view(topo);

    const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
    const double   dt        = sim.time_step_size();

    gpu::launch_leap_frog_velocity(view.old().vel, view.old().force,
                                    view.current().vel, topo_view.mass,
                                    num_atoms, dt, stream);

    // Vouch for the velocity we just wrote: no CPU round trip, and
    // don't let a subsequent configuration_view() read re-download and
    // clobber it with the (now stale) CPU value. gpu_mirror_touches()
    // == 0 keeps Algorithm_Sequence::run()'s default invalidation from
    // erasing this immediately after apply() returns.
    sim.cuda().mark_gpu_dirty(conf, gpu::MIRROR_VEL, stream);

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

    // Own stream: if a GPU-native thermostat (Berendsen/NoseHoover<
    // gpuBackend>) ran between Velocity and Position and rescaled VEL on
    // its own stream, this lets configuration_view() wait on that
    // producer's event instead of relying on legacy-default-stream
    // ordering.
    static cudaStream_t stream = 0;
    if (stream == 0) cudaStreamCreate(&stream);

    // POS and VEL are both already marked fresh on the GPU mirror
    // (Velocity's full upload covered POS, and its mark_gpu_dirty()
    // covered VEL) -- this resolves to zero extra transfers, not a
    // hand-picked "skip the sync" boolean. If some intervening
    // CPU-side algorithm (e.g. a thermostat) had run between Velocity
    // and Position, Algorithm_Sequence::run()'s default invalidation
    // after it would have cleared these bits, and this call would
    // correctly resync from CPU instead of silently reading stale GPU
    // state.
    gpu::Configuration::View view =
        sim.cuda().configuration_view(conf, gpu::MIRROR_POS | gpu::MIRROR_VEL, stream);

    gpu::launch_leap_frog_position(view.old().pos, view.current().vel,
                                    view.current().pos, num_atoms, dt, stream);

    // Stay GPU-resident -- no eager sync-back. Vouch for the position
    // we just wrote (mark_gpu_dirty(), no CPU round trip); the existing
    // generic mechanism (Algorithm_Sequence::run()'s flush_gpu_dirty()
    // before any algorithm whose gpu_mirror_touches() includes POS/VEL
    // -- the default for every CPU-only algorithm) publishes it lazily,
    // only when something genuinely needs it, same as every other
    // GPU-native writer in this codebase. gpu_mirror_touches() == 0
    // (below) keeps that same generic invalidation from immediately
    // erasing the freshness this call just marked. The one remaining
    // *unconditional* publish point is trajectory/checkpoint writing in
    // program/md.cc, which lives outside any Algorithm_Sequence::run()
    // this deferred mechanism could hook into -- see io::Out_
    // Configuration::needs_gpu_mirror_flush().
    sim.cuda().mark_gpu_dirty(conf, gpu::MIRROR_POS | gpu::MIRROR_VEL, stream);

    this->m_timer.stop();
    return 0;
}

// explicit instantiations for linker
template class algorithm::Leap_Frog_Position<util::gpuBackend>;
template class algorithm::Leap_Frog_Velocity<util::gpuBackend>;
