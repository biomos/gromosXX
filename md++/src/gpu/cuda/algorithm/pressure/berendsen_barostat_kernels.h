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
 * @file berendsen_barostat_kernels.h
 * Host-callable entry point for the GPU-native Berendsen barostat
 * position-scaling kernel. Isotropic, anisotropic, and semi-anisotropic
 * pressure coupling all reduce to the same operation on the O(num_atoms)
 * side -- scale every position by a fixed 3x3 matrix (diagonal for the
 * first two/three modes, computed once per step on the host from the
 * already CPU-resident pressure tensor -- see Berendsen_Barostat<
 * gpuBackend>::apply()). One generic matrix-vector kernel covers all
 * three; only the (tiny, O(9), always host-side) matrix itself differs
 * per mode.
 *
 * Same host/device split as leap_frog_kernels.h: this header stays
 * host-includable (no __global__ declarations), the real kernel is
 * private to the .cu.
 */

#pragma once

#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "math/gmath.h"

namespace gpu {

  /**
   * @brief pos[i] = mu * pos[i] (3x3 matrix-vector product), in place,
   * for all atoms. Mirrors Berendsen_Barostat<cpuBackend>::apply()'s
   * per-atom scaling loop exactly -- see that function for the per-mode
   * derivation of `mu` (isotropic: uniform diagonal; anisotropic/semi-
   * anisotropic: per-axis diagonal).
   */
  void launch_barostat_scale_positions(math::CuVArray::View pos,
                                        unsigned num_atoms,
                                        const FPL9_TYPE & mu,
                                        cudaStream_t stream = 0);

}
