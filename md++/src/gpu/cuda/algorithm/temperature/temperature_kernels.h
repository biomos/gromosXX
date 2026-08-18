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
 * @file temperature_kernels.h
 * Host-callable entry points for the GPU-native temperature-group
 * velocity reduction and thermostat velocity-scaling kernels
 * (Temperature_Calculation<gpuBackend>, Berendsen_Thermostat<gpuBackend>,
 * PLAN.md §10 step 13). Same host/device split as leap_frog_kernels.h /
 * remove_com_motion_kernels.h.
 *
 * Both algorithms need a temperature group's centre-of-mass velocity,
 * which is inherently a two-level (atom -> molecule) reduction -- see
 * launch_group_velocity_reduce()'s doc comment for the exact algebra
 * that lets one kernel serve both callers.
 */

#pragma once

#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "math/gmath.h"

namespace gpu {

  /**
   * @brief Per temperature group g (`group_index[i]` for atom i, static,
   * built once from topo.temperature_groups()): `sums[5g+0] = sum(mass)`,
   * `sums[5g+1..3] = sum(mass * v)`, `sums[5g+4] = sum(mass * |v|^2)`.
   * `sums` must have at least `5 * num_groups` elements.
   *
   * Called with `vel` = `conf.current().vel` to get everything needed
   * for the thermostat's feedback (`new_com_v`, `new_ekin`) and, called
   * a second time with `vel` = `conf.old().vel`, everything needed for
   * the CPU reference's "averaged" diagnostic energies too -- by
   * linearity: `com_v_avg = 0.5*(new_com_v + old_com_v)`,
   * `ekin_avg = 0.25*(sum(mass*|new_v|^2) + sum(mass*|old_v|^2))`
   * (verified against `state_properties.cc`'s `molecular_translational_
   * ekin`; both reduce to the same per-group sums this kernel already
   * computes, no separate "averaged" kernel needed).
   *
   * Direct global atomicAdd per atom into `sums[5g..5g+4]`, not a
   * per-block shared-memory bucket-then-flush -- unlike the multi-
   * energy-group tile kernel (lj_crf_tiles.cu), the group count here is
   * one per (typically rigid, few-atom) *temperature group*, i.e.
   * essentially one per solvent molecule for a real system, not a
   * handful of user-chosen energy groups. A shared-memory version needs
   * `5*num_groups` doubles of *dynamic* shared memory per block, which
   * blows past the default ~48KB limit already for a few-thousand-atom
   * real system and fails the kernel launch outright
   * (cudaErrorInvalidValue) -- found via extended_test/ubiquitin
   * (~7045 groups needs 275KB). Direct atomics have no such ceiling,
   * and since each group is only a few atoms, per-group atomic
   * contention stays low regardless of `num_groups`.
   */
  void launch_group_velocity_reduce(math::CuVArray::View vel,
                                     const float* mass,
                                     const unsigned* group_index,
                                     unsigned num_atoms,
                                     unsigned num_groups,
                                     double* sums);

  /**
   * @brief `vel(i) = scale[com_bath_of_atom[i]] * com_v[group_index[i]]
   * + scale[ir_bath_of_atom[i]] * (vel(i) - com_v[group_index[i]])`,
   * for all atoms. In-place.
   *
   * This single formula covers both of Thermostat::scale()'s CPU-side
   * cases (jointly coupled: com_bath == ir_bath: algebraically reduces
   * to `scale[b] * vel(i)` exactly; separately coupled: the general
   * case) -- verified by hand, no per-range branching needed on the GPU
   * side at all, just three static per-atom lookups.
   *
   * `com_v` must have `num_groups` entries (FPL3_TYPE, the per-group
   * centre-of-mass velocity already reduced by
   * launch_group_velocity_reduce() + a host-side divide by that group's
   * mass). `bath_scale` must have `num_baths` entries (uploaded fresh
   * every call -- changes every step, tiny).
   */
  void launch_thermostat_scale_apply(math::CuVArray::View vel,
                                      const unsigned* group_index,
                                      const unsigned* com_bath_of_atom,
                                      const unsigned* ir_bath_of_atom,
                                      const FPL3_TYPE* com_v_per_group,
                                      const double* bath_scale,
                                      unsigned num_atoms);

}
