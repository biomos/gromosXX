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

#include "gpu/cuda/math/box.h"

namespace configuration {
    class Configuration;
}
namespace gpu {

    /**
     * @brief Holds GPU-side state copies of configuration data.
     * 
     */
    struct Configuration {
        /**
         * @brief Holds current and old state
         * 
         */
        struct State {
            mutable math::CuVArray pos;
            mutable math::CuVArray vel;
            mutable math::CuVArray force;
            mutable math::CuVArray constraint_force;
            Box* box;
            float9* virial_tensor;
            float9* kinetic_energy_tensor;
            float9* pressure_tensor;
            void* tensors_block; // keep for deallocation

            State() :
                        box(nullptr),
                        virial_tensor(nullptr),
                        kinetic_energy_tensor(nullptr),
                        pressure_tensor(nullptr)
                {

                // keep box separate
                cudaMalloc(&box, sizeof(Box));
                // allocate all tensors at once
                const size_t tensors_bytes = sizeof(*virial_tensor)
                                            + sizeof(*kinetic_energy_tensor)
                                            + sizeof(*pressure_tensor);
                cudaMalloc(&tensors_block, tensors_bytes);

                char* base              = reinterpret_cast<char*>(tensors_block);
                virial_tensor           = reinterpret_cast<float9*>(base += sizeof(float9));
                kinetic_energy_tensor   = reinterpret_cast<float9*>(base += sizeof(float9));
                pressure_tensor         = reinterpret_cast<float9*>(base += sizeof(float9));
            }
            ~State() {
                cudaFree(box);
                cudaFree(tensors_block);
            }

            /**
             * @brief Export a __host__-side view to be used in kernel call
             * const pointers to non-const device data
             * 
             */
            struct StateView {
                size_t size                   = 0;
                float3* pos                   = nullptr;
                float3* vel                   = nullptr;
                float3* force                 = nullptr;
                float3* constraint_force      = nullptr;
                Box* box                      = nullptr;
                float9* virial_tensor         = nullptr;
                float9* kinetic_energy_tensor = nullptr;
                float9* pressure_tensor       = nullptr;
            };

            /**
             * @brief cache the view 
             * 
             */
            mutable StateView cached_view;
            
            /**
             * @brief Update the view with the latest data from device memory
             * 
             */
            void update_view() const {
                cached_view = StateView{
                    .size = pos.size(),
                    .pos = pos.size() ? pos.data() : nullptr,
                    .vel = vel.size() ? vel.data() : nullptr,
                    .force = force.size() ? force.data() : nullptr,
                    .constraint_force = constraint_force.size() ? constraint_force.data() : nullptr,
                    .box = box,
                    .virial_tensor = virial_tensor,
                    .kinetic_energy_tensor = kinetic_energy_tensor,
                    .pressure_tensor = pressure_tensor
                };
            }

            /**
             * @brief Get the view of the state data.
             * check for nullptr before use
             * 
             * @return View 
             */
            StateView const& view() const {
                return cached_view;
            }

            /**
             * @brief Resize the state to hold `num_atoms` atoms.
             * 
             * @param num_atoms 
             */
            void resize(size_t num_atoms) {
                pos.resize(num_atoms);
                vel.resize(num_atoms);
                force.resize(num_atoms);
                constraint_force.resize(num_atoms);
                update_view();
            }
        };

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
         * @brief if the view was created, assume mutation
         * 
         */
        // bool mutated = false;

        /**
         * @brief Update the configuration
         * 
         */
        void copy_to_device(configuration::Configuration& conf);

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
         * @brief View of the configuration
         * 
         */
        struct View {
            /**
             * @brief view of current state
             * 
             */
            State::StateView m_current;

            /**
             * @brief view of old state
             * 
             */
            State::StateView m_old;

            /**
             * @brief get current state from host, or device
             * 
             */
            HOSTDEVICE const State::StateView& current() {
                return m_current;
            }


            /**
             * @brief get old state from host, or device
             * 
             */
            HOSTDEVICE const State::StateView& old() {
                return m_old;
            }
        };

        /**
         * @brief Create view of the configuration
         * 
         * @return View 
         */
        View view() const {
            // mutated = true;
            return View{current.view(), old.view()};
        }
    };
}