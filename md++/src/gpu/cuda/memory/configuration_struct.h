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
 * @file configuration_struct.h
 * A light-weight struct for GPU holding essential Configuration data
 * The struct holds all device pointers
 */

#pragma once

#include <utility>
#include "types.h"
#include "cuvector.h"

#include "gpu/cuda/cuhostdevice.h"
#include "gpu/cuda/math/box.h"

namespace configuration {
    class Configuration;
}
namespace gpu {
    /**
     * @brief Export a __host__-side view to be used in kernel call.
     */
    struct ConfigurationStateView {
        math::CuVArray::View pos;
        math::CuVArray::View vel;
        math::CuVArray::View force;
        math::CuVArray::View constraint_force;
        Box* box;
        FPL9_TYPE* virial_tensor;
        FPL9_TYPE* kinetic_energy_tensor;
        FPL9_TYPE* pressure_tensor;

        HOSTDEVICE ConfigurationStateView() = default;

        HOSTDEVICE ConfigurationStateView(
            math::CuVArray::View p,
            math::CuVArray::View v,
            math::CuVArray::View f,
            math::CuVArray::View cf,
            Box* b = nullptr,
            FPL9_TYPE* vt = nullptr,
            FPL9_TYPE* ket = nullptr,
            FPL9_TYPE* pt = nullptr)
        : pos(p),
        vel(v),
        force(f),
        constraint_force(cf),
        box(b),
        virial_tensor(vt),
        kinetic_energy_tensor(ket),
        pressure_tensor(pt) {}
    };

    /**
     * @brief Holds exclusive pointers to current and old state
     * 
     */
    struct ConfigurationState {
        using StateView = ConfigurationStateView;
        // GPU-side arrays
        mutable math::CuVArray pos;
        mutable math::CuVArray vel;
        mutable math::CuVArray force;
        mutable math::CuVArray constraint_force;

        Box* box = nullptr;
        FPL9_TYPE* virial_tensor = nullptr;
        FPL9_TYPE* kinetic_energy_tensor = nullptr;
        FPL9_TYPE* pressure_tensor = nullptr;
        void* tensors_block = nullptr;

        ConfigurationState() {
            // Allocate Box
            cudaMalloc(&box, sizeof(Box));
            // allocate all tensors at once
            const size_t tensors_bytes = sizeof(*virial_tensor)
                                    + sizeof(*kinetic_energy_tensor)
                                    + sizeof(*pressure_tensor);
            cudaMalloc(&tensors_block, tensors_bytes);

            char* base              = reinterpret_cast<char*>(tensors_block);
            virial_tensor           = reinterpret_cast<FPL9_TYPE*>(base);
            kinetic_energy_tensor   = reinterpret_cast<FPL9_TYPE*>(base += sizeof(FPL9_TYPE));
            pressure_tensor         = reinterpret_cast<FPL9_TYPE*>(base += sizeof(FPL9_TYPE));
        }

        ~ConfigurationState() {
            if (box) cudaFree(box);
            if (tensors_block) cudaFree(tensors_block);
        }

        /**
         * @brief Construct StateView on demand
         * 
         * @return HOSTDEVICE 
         */
        __host__ StateView view() {
            return StateView{
                pos.view(),
                vel.view(),
                force.view(),
                constraint_force.view(),
                box,
                virial_tensor,
                kinetic_energy_tensor,
                pressure_tensor
            };
        }

        void resize(size_t num_atoms) {
            pos.resize(num_atoms);
            vel.resize(num_atoms);
            force.resize(num_atoms);
            constraint_force.resize(num_atoms);
        }
    };

    /**
     * @brief View of the configuration
     * 
     */
    struct ConfigurationView {
        using StateView = ConfigurationStateView;
        /**
         * @brief view of current state
         * 
         */
        StateView m_current;

        /**
         * @brief view of old state
         * 
         */
        StateView m_old;

        /**
         * @brief get current state from host, or device
         * 
         */
        HOSTDEVICE StateView& current() { return m_current; }


        /**
         * @brief get old state from host, or device
         * 
         */
        HOSTDEVICE StateView& old() { return m_old; }
    };

    /**
     * @brief Raw device pointers for passing directly to CUDA kernels.
     * Use instead of ConfigurationStateView when the kernel needs FPL3_TYPE*.
     */
    struct ConfigurationRawPtrs {
        FPL3_TYPE* pos     = nullptr;
        FPL3_TYPE* vel     = nullptr;
        FPL3_TYPE* force   = nullptr;
        FPL3_TYPE* constraint_force = nullptr;
        Box*       box     = nullptr;
    };

    /**
     * @brief Holds GPU-side state copies of configuration data.
     *
     */
    struct Configuration {
        using State = ConfigurationState;
        using View = ConfigurationView;
        using RawPtrs = ConfigurationRawPtrs;

        /**
         * @brief CudaManager's freshness bookkeeping (gpu::MirrorField
         * bits, gpu/mirror_fields.h): which fields are currently known
         * to be valid on the GPU side without a resync from CPU. Not
         * touched by this struct's own methods -- CudaManager owns the
         * read/write of this bitmask, layered on top of the dumb copy
         * routines below.
         */
        unsigned gpu_fresh_fields = 0;

        /**
         * @brief Subset of gpu_fresh_fields that the GPU holds a value
         * for which the CPU-side configuration::Configuration has NOT
         * yet seen (set by mark_gpu_dirty(), e.g. Leap_Frog_Velocity
         * <gpuBackend>'s velocity write). Distinct from
         * gpu_fresh_fields because a field can be "fresh" (matches
         * CPU, just uploaded) without being "dirty" (GPU-only,
         * diverging from CPU) -- only dirty fields need publishing to
         * CPU before a CPU-side algorithm that might read/write them
         * runs; see CudaManager::flush_gpu_dirty().
         */
        unsigned gpu_dirty_fields = 0;

        /**
         * @brief The current state
         * 
         */
        State current;

        /**
         * @brief The old state
         * 
         */
        State old;

        /**
         * @brief Copy all arrays from CPU configuration to GPU (full sync).
         */
        void copy_to_device(configuration::Configuration& conf);

        /**
         * @brief Copy only positions and velocities from CPU to GPU (per-step).
         */
        void copy_pos_vel_to_device(const configuration::Configuration& conf);

        /**
         * @brief Copy GPU forces back to the CPU configuration.
         * Calls cudaDeviceSynchronize() internally.
         */
        void copy_forces_from_device(configuration::Configuration& conf);

        /**
         * @brief Copy GPU positions and velocities (current only -- the
         * "old" state is never GPU-computed, no need to fetch it) back
         * to the CPU configuration. Calls cudaDeviceSynchronize()
         * internally. Symmetric counterpart to copy_pos_vel_to_device(),
         * for GPU-native integrators (e.g. Leap_Frog_*<gpuBackend>) that
         * compute new positions/velocities on the GPU mirror directly
         * and need to publish the result back to the CPU-authoritative
         * state once, rather than round-tripping after every kernel.
         */
        void copy_pos_vel_from_device(configuration::Configuration& conf);

        /**
         * @brief Copy GPU constraint_force (current+old) and virial_tensor
         * (current+old) back to the CPU configuration. Calls
         * cudaDeviceSynchronize() internally. Counterpart to
         * copy_pos_vel_from_device() for the constraint algorithms
         * (SHAKE/M-SHAKE/LINCS/SETTLE), which accumulate into the
         * shared GPU-resident constraint_force/virial_tensor buffers
         * (gpu::MIRROR_CONSTRAINT_FORCE/MIRROR_VIRIAL) instead of each
         * keeping a private one and round-tripping through the CPU
         * every apply() call.
         */
        void copy_constraint_data_from_device(configuration::Configuration& conf);

        /**
         * @brief Allow to exchange states efficiently
         * Views are exchanged implicitly as well
         * 
         */
        void exchange_state() {
            std::swap(current, old);
        }

        /**
         * @brief resize the arrays
         * 
         * @param num_atoms 
         */
        void resize(size_t num_atoms) {
            current.resize(num_atoms);
            old.resize(num_atoms);
        }

        /**
         * @brief Create view of the configuration
         *
         * @return View
         */
        View view() { return View{current.view(), old.view()}; }

        /**
         * @brief Raw device pointers of the current state, for direct kernel use.
         */
        RawPtrs current_raw() {
            return { current.pos.data(), current.vel.data(),
                     current.force.data(), current.constraint_force.data(),
                     current.box };
        }
    };
}