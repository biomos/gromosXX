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
 * @file periodicity.h
 * implementation of the periodic boundary condition functions for GPU.
 */

#pragma once

#undef MODULE
#undef SUBMODULE
#define MODULE gpu
#define SUBMODULE math

#include "cuhostdevice.h"
#include "math/box.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/math/cumath.h"
// #include "gpu/cuda/memory/topology_struct.h"
// #include "gpu/cuda/memory/configuration_struct.h"
// #include "math/periodicity.h"
// #include "gpu/cuda/kernels/hello_world.h"

#include "math/morton.h"

namespace configuration {
  class Configuration;
}

namespace topology {
  class Topology;
}

namespace gpu {
/**
 * @class Periodicity
 * GPU-friendly periodic boundary condition functions.
 * Conceptually mirrors math::Periodicity, but designed for device.
 */
    template <math::boundary_enum BOUNDARY>
    class Periodicity {
        public:
            // Constructor initializes device pointers
            __host__ Periodicity(math::Box box) {
                m_box = gpu::Box(box);
                m_half_box.x = 0.5 * abs(m_box(0));
                m_half_box.y = 0.5 * abs(m_box(1));
                m_half_box.z = 0.5 * abs(m_box(2));
            }

            // Device-friendly function to put vector into box
            template <typename VecType>
            HOSTDEVICE void put_into_box(VecType &v) const {
                VecType o{0., 0., 0.};
                v = this->nearest_image(v, o);
            }

            /**
             * @brief Computes nearest image of vector v1 with respect to v2
             * 
             * @tparam VecType 
             * @param v1 
             * @param v2 
             * @return HOSTDEVICE 
             */
            template <typename VecType>
            HOSTDEVICE VecType nearest_image(const VecType &v1,
                                            const VecType &v2) const {
                VecType nim;
                
                if constexpr (BOUNDARY == math::vacuum) {
                nim = v1 - v2;
                } else if constexpr (BOUNDARY == math::rectangular) {
                nim.x = v1.x - v2.x;
                while (nim.x > m_half_box.x) nim.x -= 2.0 * m_half_box.x;
                while (nim.x < -m_half_box.x) nim.x += 2.0 * m_half_box.x;

                nim.y = v1.y - v2.y;
                while (nim.y > m_half_box.y) nim.y -= 2.0 * m_half_box.y;
                while (nim.y < -m_half_box.y) nim.y += 2.0 * m_half_box.y;

                nim.z = v1.z - v2.z;
                while (nim.z > m_half_box.z) nim.z -= 2.0 * m_half_box.z;
                while (nim.z < -m_half_box.z) nim.z += 2.0 * m_half_box.z;

                } else {
                nim.x = nim.y = nim.z = 0; // avoid compiler warnings
                assert(false);
                }
                return nim;
            }

            /**
             * @brief Set the cell size based on the cutoff and box size
             * We make sure the cell size > cutoff, so the halo always
             * contains only the adjacent cells
             * 
             * @param box 
             * @param max_cutoff 
             */
            __host__ void set_cell_size(double max_cutoff) {
                if constexpr (BOUNDARY == math::vacuum) {
                    m_num_cells = {1,1,1};
                    m_inv_cell_size = {0,0,0};
                    
                } else if constexpr (BOUNDARY == math::rectangular) {
                    double3 box{abs(m_box(0)),
                                abs(m_box(1)),
                                abs(m_box(2))};
                    double imax_cutoff = 1 / max_cutoff;
                    m_num_cells.x = (unsigned)floor(box.x * imax_cutoff);
                    m_num_cells.y = (unsigned)floor(box.y * imax_cutoff);
                    m_num_cells.z = (unsigned)floor(box.z * imax_cutoff);

                    // set to >=1
                    m_num_cells.x = max(1, m_num_cells.x);
                    m_num_cells.y = max(1, m_num_cells.y);
                    m_num_cells.z = max(1, m_num_cells.z);

                    // calculate inverted sizes, will be used to assign entities to cells
                    m_inv_cell_size.x = m_num_cells.x / box.x;
                    m_inv_cell_size.y = m_num_cells.y / box.y;
                    m_inv_cell_size.z = m_num_cells.z / box.z;
                    
                } else {
                    std::cerr << "This boundary type is not implemented in gpu::Periodicity" << std::endl;
                }
                return;
            }

            /**
             * @brief get cell indices,
             * call this only after using put_into_box(pos)!
             * 
             * @param pos 
             * @return HOSTDEVICE 
             */
            HOSTDEVICE ushort4 get_cell(FPL3_TYPE& pos) const {
                // compute indices
                unsigned ix = (unsigned)floorf((pos.x + m_half_box.x) * m_inv_cell_size.x);
                unsigned iy = (unsigned)floorf((pos.y + m_half_box.y) * m_inv_cell_size.y);
                unsigned iz = (unsigned)floorf((pos.z + m_half_box.z) * m_inv_cell_size.z);
                
                // wrap edge cases
                ix = min(ix, m_num_cells.x - 1);
                iy = min(iy, m_num_cells.y - 1);
                iz = min(iz, m_num_cells.z - 1);

                assert(
                ix < 32 &&
                iy < 32 &&
                iz < 32 ); // we can fit only 3x5 bits -> max int is 32
                
                ushort4 cell;
                cell.x = (ushort)ix;
                cell.y = (ushort)iy;
                cell.z = (ushort)iz;
                // also generate the morton index in case we want to have more efficient ordering
                cell.w = (ushort)morton_index_16(ix, iy, iz);

                return cell;
            }

            /**
             * @brief Prepare chargegroups, i.e. calculate centers of geometry (cog), and translate atoms into box
             * optional:
             *      - store cog in provided array
             *      - store cell id in provided array
             * 
             * @tparam ArrT array of positions, math::CuVArray
             * @tparam CGCogArr pointer to float3 or double3 array
             * @tparam CGCellArr pointer to ushort4 array
             * @param cg_i 
             * @param cg_begin 
             * @param cg_end 
             * @param pos 
             * @param num_solute_chargegroups 
             * @param cg_cog_ptr 
             * @param cg_cells_ptr 
             * @return __device__ 
             */
            template <typename ArrT, typename CGCogArr = void, typename CGCellArr = void>
            __device__ void prepare_chargegroup(
                unsigned cg_i,
                unsigned cg_begin,
                unsigned cg_end,
                ArrT& pos,
                unsigned num_solute_chargegroups = 0,
                CGCogArr* cg_cog_ptr = nullptr,
                CGCellArr* cg_cells_ptr = nullptr
            ) const {
                // compute cog
                FPL3_TYPE cog = (cg_i < num_solute_chargegroups)
                                ? calculate_centre_of_geometry(pos, cg_begin, cg_end)
                                : pos(cg_begin);

                // put into box
                FPL3_TYPE v_box = cog;
                put_into_box(v_box);

                // optional: store cog
                if (cg_cog_ptr && cg_i < num_solute_chargegroups) {
                    (*cg_cog_ptr)(cg_i) = v_box;
                }

                // optional: store cell
                if (cg_cells_ptr) {
                    (*cg_cells_ptr)(cg_i) = get_cell(v_box);
                }

                // translate positions
                FPL3_TYPE trans = v_box - cog;
                for (unsigned a_i = cg_begin; a_i < cg_end; ++a_i) {
                    pos(a_i) += trans;
                }
            }

        private:
            /**
             * @brief full box struct
             * 
             */
            gpu::Box m_box;

            /**
             * @brief half box
             * 
             */
            FPL3_TYPE m_half_box;

            /**
             * @brief inverted cell size
             * 
             */
            FPL3_TYPE m_inv_cell_size;

            /**
             * @brief number of cells in each dimension
             * 
             */
            dim3 m_num_cells;
    };

} // namespace gpu