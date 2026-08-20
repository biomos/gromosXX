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
 * @file remove_com_motion_kernels.h
 * Host-callable entry points for the GPU-native centre-of-mass motion
 * removal kernels (Remove_COM_Motion<gpuBackend>, PLAN.md §10 step 13 --
 * broadening algorithm coverage). Same host/device split as
 * leap_frog_kernels.h: this header stays host-includable (no
 * __global__ declarations), the real kernels are private to the .cu.
 *
 * All reductions produce `double` sums (not FPL_TYPE) regardless of
 * build precision, same rationale as the nonbonded energy/virial
 * reductions: summing thousands of per-atom contributions into a
 * single accumulator loses less precision in double, and the CPU
 * reference (remove_com_motion_cpu.cc) itself accumulates in double
 * throughout -- matching that, not FPL_TYPE, is what makes the
 * CPU/GPU comparison test meaningful.
 */

#pragma once

#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "math/gmath.h"

namespace gpu {

  /**
   * @brief sums[0..2] = sum_i(mass[i] * vel[i]) (total momentum),
   * sums[3] = sum_i(mass[i]) (total mass). `sums` must have at least 4
   * elements, caller-owned (reused across calls, not allocated here).
   * Does not synchronize -- caller must cudaDeviceSynchronize() (or
   * check a stream/event) before reading `sums`.
   */
  void launch_com_translation_reduce(math::CuVArray::View vel,
                                      const float* mass,
                                      unsigned num_atoms,
                                      double* sums,
                                      cudaStream_t stream = 0);

  /**
   * @brief vel[i] -= com_v, for all atoms. In-place.
   */
  void launch_com_translation_apply(math::CuVArray::View vel,
                                     double com_v_x, double com_v_y, double com_v_z,
                                     unsigned num_atoms,
                                     cudaStream_t stream = 0);

  /**
   * @brief sums[0..2] = sum_i(mass[i] * vel[i]), sums[3..5] =
   * sum_i(mass[i]*pos[i] - 0.5*mass[i]*vel[i]*dt), sums[6] =
   * sum_i(mass[i]). `sums` must have at least 7 elements. First pass
   * of Remove_COM_Motion<gpuBackend>::remove_com_rotation() -- produces
   * com_v/com_r (after the host divides by total mass), needed before
   * the second pass (launch_com_rotation_reduce_pass2) can run.
   */
  void launch_com_rotation_reduce_pass1(math::CuVArray::View pos,
                                         math::CuVArray::View vel,
                                         const float* mass,
                                         double dt,
                                         unsigned num_atoms,
                                         double* sums,
                                         cudaStream_t stream = 0);

  /**
   * @brief sums[0..2] = angular momentum L = sum_i(mass[i] *
   * cross(r_i, vel[i] - com_v)), sums[3..11] = inertia tensor I
   * (flattened row-major, index i*3+j matching math::Matrix::
   * operator()(i,j)), where r_i = pos[i] - 0.5*dt*vel[i] - com_r.
   * `sums` must have at least 12 elements. Exact CPU formula
   * (remove_com_motion_cpu.cc's remove_com_rotation).
   */
  void launch_com_rotation_reduce_pass2(math::CuVArray::View pos,
                                         math::CuVArray::View vel,
                                         const float* mass,
                                         double dt,
                                         double com_v_x, double com_v_y, double com_v_z,
                                         double com_r_x, double com_r_y, double com_r_z,
                                         unsigned num_atoms,
                                         double* sums,
                                         cudaStream_t stream = 0);

  /**
   * @brief vel[i] -= cross(com_O, r_i), r_i = pos[i] - 0.5*dt*vel[i] -
   * com_r (recomputed here, not stored, from the same pos/vel/dt/com_r
   * already used by pass2 -- avoids a separate per-atom r_i buffer).
   * In-place.
   */
  void launch_com_rotation_apply(math::CuVArray::View pos,
                                  math::CuVArray::View vel,
                                  double dt,
                                  double com_r_x, double com_r_y, double com_r_z,
                                  double com_O_x, double com_O_y, double com_O_z,
                                  unsigned num_atoms,
                                  cudaStream_t stream = 0);

}
