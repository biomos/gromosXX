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
 * @file leap_frog.h
 * Leap-frog integration algorithms.
 *
 * Both Leap_Frog_Velocity and Leap_Frog_Position are templated on Backend:
 *   - util::cpuBackend : existing OpenMP-parallelised CPU implementation
 *   - util::gpuBackend : CUDA kernel implementation (operates on GPU arrays)
 *
 * Usage with backend dispatch:
 * @code
 *   auto* lfv = algorithm::make_algorithm<Leap_Frog_Velocity>(sim);
 *   auto* lfp = algorithm::make_algorithm<Leap_Frog_Position>(sim);
 * @endcode
 */

#pragma once

namespace algorithm
{

// ─────────────────────────────────────────────────────────────────────────────

/**
 * @class Leap_Frog_Velocity
 * v(t+dt/2) = v(t-dt/2) + dt/m * F(t)
 */
template <typename Backend = util::cpuBackend>
class Leap_Frog_Velocity : public Algorithm, private AlgorithmB<Backend>
{
public:
    template <typename B>
    static constexpr bool is_supported_backend =
        std::is_same_v<B, util::cpuBackend> ||
        std::is_same_v<B, util::gpuBackend>;

    Leap_Frog_Velocity() : Algorithm("Leap_Frog_Velocity") {}
    virtual ~Leap_Frog_Velocity() {}

    virtual int apply(topology::Topology& topo,
                      configuration::Configuration& conf,
                      simulation::Simulation& sim);

    // gpuBackend: already tracks the fields it touches precisely via
    // sim.cuda() (mark_gpu_dirty() for VEL, a full resync for
    // everything else) -- returning 0 stops Algorithm_Sequence::run()'s
    // default post-apply() invalidation from immediately erasing the
    // VEL freshness this just set, which Leap_Frog_Position<gpuBackend>
    // depends on reading without a resync. cpuBackend: matches the
    // base class default, harmless (nothing GPU-resident to protect).
    virtual unsigned gpu_mirror_touches() const override {
        if constexpr (std::is_same_v<Backend, util::gpuBackend>)
            return 0u;
        else
            return gpu::MIRROR_ALL;
    }

    virtual int init(topology::Topology& topo,
                     configuration::Configuration& conf,
                     simulation::Simulation& sim,
                     std::ostream& os    = std::cout,
                     bool          quiet = false)
    {
        if (!quiet) {
            if constexpr (std::is_same_v<Backend, util::gpuBackend>)
                os << "INTEGRATION\n\tLeap frog velocity (GPU)\n";
            else
                os << "INTEGRATION\n\tLeap frog velocity\n";
        }
        return 0;
    }
};

// ─────────────────────────────────────────────────────────────────────────────

/**
 * @class Leap_Frog_Position
 * x(t+dt) = x(t) + dt * v(t+dt/2)
 */
template <typename Backend = util::cpuBackend>
class Leap_Frog_Position : public Algorithm, private AlgorithmB<Backend>
{
public:
    template <typename B>
    static constexpr bool is_supported_backend =
        std::is_same_v<B, util::cpuBackend> ||
        std::is_same_v<B, util::gpuBackend>;

    Leap_Frog_Position() : Algorithm("Leap_Frog_Position") {}
    virtual ~Leap_Frog_Position() {}

    virtual int apply(topology::Topology& topo,
                      configuration::Configuration& conf,
                      simulation::Simulation& sim);

    // gpuBackend: never writes configuration::Configuration directly
    // (reads/writes only through sim.cuda(), and its own
    // sync_configuration_from_device() call already publishes the
    // final result to CPU before returning) -- returning 0 avoids a
    // pointless flush-before-self on its own before-apply() hook, and
    // the subsequent invalidate-after is harmless since nothing is
    // left dirty by the time it returns. cpuBackend: base class
    // default, harmless.
    virtual unsigned gpu_mirror_touches() const override {
        if constexpr (std::is_same_v<Backend, util::gpuBackend>)
            return 0u;
        else
            return gpu::MIRROR_ALL;
    }

    virtual int init(topology::Topology& topo,
                     configuration::Configuration& conf,
                     simulation::Simulation& sim,
                     std::ostream& os    = std::cout,
                     bool          quiet = false)
    {
        if (!quiet) {
            if constexpr (std::is_same_v<Backend, util::gpuBackend>)
                os << "\tLeap frog position (GPU)\nEND\n";
            else
                os << "\tLeap frog position\nEND\n";
        }
        return 0;
    }
};

} // namespace algorithm
