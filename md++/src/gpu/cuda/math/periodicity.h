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
// #include "gpu/cuda/memory/topology_struct.h"
// #include "gpu/cuda/memory/configuration_struct.h"
// #include "math/periodicity.h"
// #include "gpu/cuda/kernels/hello_world.h"

namespace configuration {
  class Configuration;
}

namespace topology {
  class Topology;
}

namespace gpu {
  template <typename VecType>
  __device__ VecType calculate_cog(VecType* pos, int begin, int end);

  /**
   * @class Periodicity
   * GPU-friendly periodic boundary condition functions.
   * Conceptually mirrors math::Periodicity, but designed for device.
   */
  template <math::boundary_enum BOUNDARY>
  struct Periodicity {
      /**
       * @brief gpu::Box accessible from host, passable as argument to device
       * 
       */
      gpu::Box m_box;

      /**
       * @brief half box
       * 
       */
      FPL3_TYPE m_half_box;

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

      // Optional: nearest image calculation
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
            nim.x = 0;
            nim.y = 0;
            nim.z = 0;
            assert(false);
          }
          return nim;
      }

      // Additional GPU helper functions can go here (translations, shifts)
  };

} // namespace gpu