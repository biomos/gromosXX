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
        std::is_same_v<B, util::cpuBackend>;
    // TODO(cleanup): gpuBackend disabled for now -- integration/leap_frog_gpu.cu
    // exists but is not wired into any CMakeLists.txt and calls GPU-mirror APIs
    // (conf.m_gpu->old_raw(), conf.copy_pos_vel_to_gpu()) that don't fully exist
    // (gpu::Configuration has no old_raw()) and are being removed anyway per
    // PLAN.md §3.2. Claiming gpuBackend support here without a linkable
    // definition breaks the cuda-on link. Re-enable once PLAN.md §10 roadmap
    // step 11 re-ports leap_frog_gpu.cu against sim.cuda().configuration_view().
    // ||  std::is_same_v<B, util::gpuBackend>;

    Leap_Frog_Velocity() : Algorithm("Leap_Frog_Velocity") {}
    virtual ~Leap_Frog_Velocity() {}

    virtual int apply(topology::Topology& topo,
                      configuration::Configuration& conf,
                      simulation::Simulation& sim);

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
        std::is_same_v<B, util::cpuBackend>;
    // TODO(cleanup): see Leap_Frog_Velocity above -- same reason, same plan.
    // ||  std::is_same_v<B, util::gpuBackend>;

    Leap_Frog_Position() : Algorithm("Leap_Frog_Position") {}
    virtual ~Leap_Frog_Position() {}

    virtual int apply(topology::Topology& topo,
                      configuration::Configuration& conf,
                      simulation::Simulation& sim);

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
