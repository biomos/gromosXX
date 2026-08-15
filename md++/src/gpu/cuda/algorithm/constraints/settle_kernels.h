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
 * @file settle_kernels.h
 * Host-callable entry point for the SETTLE kernel (PLAN.md §10 step
 * 18, constraints) -- the analytical (closed-form, non-iterative)
 * rigid 3-site-water constraint algorithm (Miyamoto & Kollman, J.
 * Comput. Chem 13, 1992), matching algorithm::Settle::solvent()
 * (settle.cc) line-for-line.
 *
 * Unlike SHAKE, SETTLE has no iteration at all -- one thread per
 * molecule does exactly the same fixed sequence of closed-form vector
 * algebra as the CPU reference, so (unlike solute SHAKE's Jacobi
 * scheme) this port *is* bit-comparable to the CPU per molecule, up to
 * ordinary floating-point non-associativity. Molecules are fully
 * independent (each owns 3 disjoint atoms), so this is embarrassingly
 * parallel across molecules -- the same shape as solvent SHAKE
 * (shake_kernels.h), just non-iterative.
 *
 * Deliberately uses plain `double`/`double3` throughout, not
 * `FPL_TYPE`/`FPL3_TYPE` -- same rationale as shake_kernels.h: the
 * position correction here directly becomes the position the entire
 * next MD step integrates from.
 *
 * v1 scope matches algorithm::Settle::init()'s own gates (checked
 * host-side in CUDA_Settle::init(), not repeated here): exactly one
 * solvent type, exactly 3 atoms per molecule, H1/H2 same mass, exactly
 * 3 distance constraints with the two O-H constraints sharing one
 * length.
 */

#pragma once

#include "gpu/cuda/memory/cuvector.h"
#include "math/gmath.h"

namespace gpu {

  /**
   * @param pos current positions, flattened double3, device, updated
   *   in place for the solvent range only.
   * @param old_pos reference (previous-step) positions, same layout,
   *   read-only.
   * @param vel current velocities, updated in place if `do_velocity`
   *   (ignored otherwise -- may be null).
   * @param first_atom global index of molecule 0's first atom
   *   (== topo.num_solute_atoms()).
   * @param num_molecules number of (single) solvent-type molecules.
   * @param mass_O, mass_H, dist_OH, dist_HH: exactly
   *   algorithm::Settle::solvent()'s inputs of the same name.
   * @param dt_i 1 / sim.time_step_size().
   * @param do_velocity whether to update `vel` (matches the CPU's
   *   `!stochastic.sd && !minimise.ntem && !analyze.analyze` gate,
   *   evaluated host-side).
   * @param constraint_force written into (device, global, solvent-
   *   range-sized; disjoint per molecule, no atomics needed) -- NOT
   *   accumulated (`=`, matching the CPU's `cons_force[k] = ...`, no
   *   pre-existing value to add to).
   * @param virial 9-double device global accumulator, atomicAdd'd
   *   (contended across molecules) -- always computed; whether it's
   *   added into `conf.old().virial_tensor` is the host's decision
   *   (matches the CPU's live `pcouple.virial == atomic_virial` gate).
   * @param error_flag device int, atomicExch'd to 1 if any molecule's
   *   geometry is degenerate (matches the CPU's three `sin(x) > 1.0`
   *   checks, which currently error but skip only that molecule via
   *   `continue`) -- caller zeroes before launch and checks after.
   */
  void launch_settle(
      double3* pos,
      const double3* old_pos,
      double3* vel,
      unsigned first_atom,
      unsigned num_molecules,
      double mass_O,
      double mass_H,
      double dist_OH,
      double dist_HH,
      double dt_i,
      bool do_velocity,
      double3* constraint_force,
      double* virial,
      int* error_flag,
      cudaStream_t stream = 0);

} // namespace gpu
