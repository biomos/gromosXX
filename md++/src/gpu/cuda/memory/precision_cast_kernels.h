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
 * @file precision_cast_kernels.h
 * Device-to-device FPL3_TYPE <-> double3 element-wise cast, over a
 * contiguous [first, first + count) atom range. For GPU-native
 * algorithms that deliberately keep their own math in full double
 * precision (e.g. CUDA_Settle -- see settle_kernels.h's doc comment)
 * while still reading/writing the shared GPU mirror, which stores
 * pos/vel at FPL precision (float under FP_PRECISION 1/2). Runs
 * entirely on-device -- no CPU round trip, unlike converting via a
 * host-side upload/download.
 */

#pragma once

#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/memory/precision.h"
#include "math/gmath.h"

namespace gpu {

  void launch_cast_fpl3_to_double3(
      const FPL3_TYPE* src, double3* dst, unsigned first, unsigned count,
      cudaStream_t stream = 0);

  void launch_cast_double3_to_fpl3(
      const double3* src, FPL3_TYPE* dst, unsigned first, unsigned count,
      cudaStream_t stream = 0);

} // namespace gpu
