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
 * @file leap_frog_kernels.h
 * Host-callable entry points for the GPU-native leap-frog integration
 * kernels (PLAN.md §10 roadmap step 11). No boundary/periodicity
 * handling is needed here -- unlike the pairlist kernels, these are
 * plain per-atom elementwise updates; box rewrapping happens elsewhere
 * (chargegroup preparation in the nonbonded path), not in integration.
 *
 * Same host/device split as lj_crf_tiles.h / displacement.h: this
 * header stays host-includable (no __global__ declarations), the real
 * kernels are private to the .cu.
 */

#pragma once

#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "math/gmath.h"

namespace gpu {

  /**
   * @brief v_new[i] = v_old[i] + f_old[i] * dt / mass[i], for all atoms.
   * Mirrors Leap_Frog_Velocity<cpuBackend>::apply()'s per-atom loop
   * exactly (leap_frog.cc). Caller is responsible for having already
   * performed the CPU-side conf.exchange_state() and having resynced
   * the GPU mirror (full resync -- old().force is needed fresh) before
   * calling this.
   */
  void launch_leap_frog_velocity(math::CuVArray::View old_vel,
                                  math::CuVArrayH::View old_force,
                                  math::CuVArray::View new_vel,
                                  const float* mass,
                                  unsigned num_atoms,
                                  double dt,
                                  cudaStream_t stream = 0);

  /**
   * @brief x_new[i] = x_old[i] + v_current[i] * dt, for all atoms.
   * Mirrors Leap_Frog_Position<cpuBackend>::apply()'s per-atom loop
   * exactly. `v_current` is read from the same GPU mirror that
   * launch_leap_frog_velocity() just wrote into -- no CPU round trip
   * between the two algorithms.
   */
  void launch_leap_frog_position(math::CuVArray::View old_pos,
                                  math::CuVArray::View current_vel,
                                  math::CuVArray::View new_pos,
                                  unsigned num_atoms,
                                  double dt,
                                  cudaStream_t stream = 0);

}
