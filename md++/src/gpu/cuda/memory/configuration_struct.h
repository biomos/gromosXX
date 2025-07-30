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
 */

#include <utility>
#include "types.h"
#include "cuvector.h"

namespace gpu {

  /**
   * @brief Holds GPU-side state copies of configuration data.
   * 
   */
    struct Configuration {
        struct State {
            math::CuVArray pos;
            math::CuVArray vel;
            math::CuVArray force;
            math::CuVArray constraint_force;
            float9* box;
            float9* virial_tensor;
            float9* kinetic_energy_tensor;
            float9* pressure_tensor;

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
            void update_view() {
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

            View const& view() const {
                return cached_view;
            }

            void resize(size_t num_atoms) {
                pos.resize(num_atoms);
                vel.resize(num_atoms);
                force.resize(num_atoms);
                constraint_force.resize(num_atoms);
                update_view();
            }
        };
        typedef State::View StateView;

        State current;
        State old;

        void exchange_state() {
            std::swap(current, old);
        }

        void resize(size_t num_atoms) {
            current.resize(num_atoms);
            old.resize(num_atoms);
        }

        // This struct is safe to copy and use on device
        class View {
            private:
                StateView m_current;
                StateView m_old;
            public:
                View(StateView current, StateView old) :
                        m_current(current), m_old(old) {}
                __host__ __device__ StateView current() const {
                    return m_current;
                }
                __host__ __device__ StateView old() const {
                    return m_old;
                }
        };

        View view() const {
            return View(current.view(), old.view());
        }
    };
}