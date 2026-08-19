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
 * @file vec3_convert.h
 * Bulk host<->device conversion between `math::Vec` (`GenericVec<double>`,
 * a plain `double d_v[3]`) and CUDA's `double3` -- used by the
 * double-precision constraint algorithms (SHAKE/SETTLE/LINCS,
 * gpu/cuda/algorithm/constraints/*) that upload/download the *entire*
 * position array every call, unlike everything else in the GPU
 * pipeline (which uses the float `FPL3_TYPE` mirror -- see PLAN.md's
 * precision note).
 *
 * `math::Vec` and `double3` are both exactly three tightly-packed
 * `double`s with identical alignment, so the whole array can be copied
 * with one `memcpy` instead of a host loop that reconstructs each
 * element field-by-field -- the loop was doing three bounds-checked
 * accessor calls (`pos(i)(0)`, `(1)`, `(2)`) and a struct construction
 * per atom, all to move bytes that are already laid out identically.
 * `static_assert`ed below so a future change to either type's layout
 * fails to compile instead of silently corrupting positions.
 */

#pragma once

#include <cstring>
#include <cuda_runtime.h>

#include "math/gmath.h"
#include "gpu/cuda/memory/precision.h"

namespace gpu {

  static_assert(sizeof(double3) == sizeof(math::Vec),
                "double3 and math::Vec must have identical layout for "
                "the bulk memcpy in vec3_convert.h to be valid");
  static_assert(alignof(double3) == alignof(math::Vec),
                "double3 and math::Vec must have identical alignment for "
                "the bulk memcpy in vec3_convert.h to be valid");

  /** Copies `n` consecutive `math::Vec`s starting at `*src` into `dst`. */
  inline void vec3_upload(double3* dst, const math::Vec* src, std::size_t n) {
    std::memcpy(dst, src, n * sizeof(double3));
  }

  /** Copies `n` consecutive `double3`s starting at `*src` into `dst`. */
  inline void vec3_download(math::Vec* dst, const double3* src, std::size_t n) {
    std::memcpy(dst, src, n * sizeof(double3));
  }

  /**
   * FPL3_TYPE counterparts for the mixed-precision constraint kernels
   * (SHAKE/M-SHAKE/SETTLE) -- unlike the double3 versions above, this
   * can't be a bulk memcpy: FPL3_TYPE is `float3` under FP_PRECISION
   * 1/2, a different byte layout than `math::Vec`'s three `double`s, so
   * each component is narrowed/widened element-by-element. Only used
   * at each apply() call's upload/download boundary, not per-iteration.
   */
  inline void vec3_upload_fpl(FPL3_TYPE* dst, const math::Vec* src, std::size_t n) {
    for (std::size_t i = 0; i < n; ++i) {
      dst[i] = make_FPL3(static_cast<FPL_TYPE>(src[i](0)),
                          static_cast<FPL_TYPE>(src[i](1)),
                          static_cast<FPL_TYPE>(src[i](2)));
    }
  }

  inline void vec3_download_fpl(math::Vec* dst, const FPL3_TYPE* src, std::size_t n) {
    for (std::size_t i = 0; i < n; ++i) {
      dst[i] = math::Vec(static_cast<double>(src[i].x),
                          static_cast<double>(src[i].y),
                          static_cast<double>(src[i].z));
    }
  }

} // namespace gpu
