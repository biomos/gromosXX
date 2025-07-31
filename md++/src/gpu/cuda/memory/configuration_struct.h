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
            float9* box;
            float9* virial_tensor;
            float9* kinetic_energy_tensor;
            float9* pressure_tensor;
            void* memory_block; // keep for deallcation
            
            /**
             * @brief size to allocate using cudaMalloc
             * 
             */
            static constexpr size_t size_in_bytes =
                4 * sizeof(float9);  // box, virial_tensor, kinetic_energy_tensor, pressure_tensor

            State() :
                        box(nullptr),
                        virial_tensor(nullptr),
                        kinetic_energy_tensor(nullptr),
                        pressure_tensor(nullptr)
                {

                // keep box separate
                cudaMalloc(&box, sizeof(float9));
                cudaMalloc(&memory_block, sizeof(float9)*3);

                char* base              = reinterpret_cast<char*>(memory_block);
                virial_tensor           = reinterpret_cast<float9*>(base += sizeof(float9));
                kinetic_energy_tensor   = reinterpret_cast<float9*>(base += sizeof(float9));
                pressure_tensor         = reinterpret_cast<float9*>(base += sizeof(float9));
            }
            ~State() {
                cudaFree(box);
                cudaFree(memory_block);
            }

            /**
             * @brief Export a __host__-side view to be used in kernel call
             * const pointers to non-const device data
             * 
             */
            struct View {
                float3* pos;
                float3* vel;
                float3* force;
                float3* constraint_force;
                float9* box;
                float9* virial_tensor;
                float9* kinetic_energy_tensor;
                float9* pressure_tensor;
            };

            /**
             * @brief cache the view 
             * 
             */
            mutable View cached_view;
            
            /**
             * @brief Update the view with the latest data from device memory
             * 
             */
            void update_view() const {
                cached_view = View{
                    pos.size() ? pos.data() : nullptr,
                    vel.size() ? vel.data() : nullptr,
                    force.size() ? force.data() : nullptr,
                    constraint_force.size() ? constraint_force.data() : nullptr,
                    box,
                    virial_tensor,
                    kinetic_energy_tensor,
                    pressure_tensor
                };
            }

            /**
             * @brief Get the view of the state data.
             * 
             * @return View 
             */
            View const& view() const {
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
            State::View m_current;

            /**
             * @brief view of old state
             * 
             */
            State::View m_old;

            /**
             * @brief get current state from host, or device
             * 
             */
            __host__ __device__ const State::View& current() {
                return m_current;
            }


            /**
             * @brief get old state from host, or device
             * 
             */
            __host__ __device__ const State::View& old() {
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