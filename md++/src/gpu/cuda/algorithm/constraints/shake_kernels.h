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
 * @file shake_kernels.h
 * Host-callable entry point for the per-solvent-molecule SHAKE kernel
 * (PLAN.md §10 step 15, constraints). One CUDA thread per solvent
 * molecule, doing the *exact same* sequential Gauss-Seidel iteration
 * (with the skip_now/skip_next per-atom convergence-tracking
 * optimization) as the CPU reference (algorithm::Shake::shake_iteration
 * / algorithm::Shake::solvent, shake.h) -- molecules are independent of
 * each other (no shared atoms), so parallelizing *across* molecules
 * while keeping each molecule's own inner loop sequential (exactly like
 * the CPU) is both correct and embarrassingly parallel. This is why
 * solvent SHAKE, unlike solute SHAKE, ports faithfully: solute SHAKE's
 * CPU reference is one large in-place Gauss-Seidel sweep over every
 * solute distance constraint with genuine cross-constraint data
 * dependencies within a single iteration -- not attempted here (see
 * cuda_shake.h's doc comment for the v1 scope gate).
 *
 * Deliberately uses plain `double`/`double3` throughout, not
 * `FPL_TYPE`/`FPL3_TYPE` -- unlike a force/energy kernel, SHAKE's
 * position corrections directly become the positions the *entire next
 * MD step* integrates from, so precision loss here is far more
 * consequential than in an energy tolerance comparison (PLAN.md's
 * quartic-bond precision note explains the magnitude for a similarly
 * cancellation-prone formula: `diff = r0^2 - dist2`). Positions are
 * uploaded/downloaded as raw doubles each call, independent of the
 * `FP_PRECISION` build setting and the float position mirror used
 * elsewhere in the GPU pipeline.
 */

#pragma once

#include "gpu/cuda/memory/cuvector.h"
#include "math/gmath.h"

namespace gpu {

  /**
   * Per-molecule skip_now/skip_next arrays are thread-local fixed-size
   * arrays sized to this constant -- solvent molecules (water, methanol,
   * etc.) are always small. Checked against the real
   * topo.solvent(i).num_atoms() at CUDA_Shake::init() time; a solvent
   * type exceeding this hard-errors rather than silently truncating.
   */
  constexpr unsigned MAX_SHAKE_ATOMS_PER_MOLECULE = 8;

  /**
   * One entry per distance constraint within a solvent molecule,
   * local (0-indexed within the molecule) atom indices.
   */
  struct ShakeConstraint {
    unsigned i;
    unsigned j;
    double r0sq;
  };

  /**
   * @param pos current positions, flattened double3, device, updated
   *   in place for solvent atoms only (indices [first_atom, first_atom
   *   + num_molecules*num_atoms_per_molecule)).
   * @param old_pos reference (previous-step) positions, same layout,
   *   read-only.
   * @param constraints per-molecule-type constraint list (same for
   *   every molecule of this solvent type), device.
   * @param num_constraints length of `constraints`.
   * @param inv_mass_local inverse mass per local atom index (device,
   *   length num_atoms_per_molecule).
   * @param num_atoms_per_molecule atoms in one molecule of this
   *   solvent type -- must be <= gpu::MAX_SHAKE_ATOMS_PER_MOLECULE
   *   (checked at host init() time).
   * @param first_atom global index of molecule 0's first atom.
   * @param num_molecules number of molecules of this solvent type.
   * @param constraint_force accumulated into (device, global,
   *   num_atoms-sized, NOT zeroed by this kernel -- caller zeroes
   *   once per call). Disjoint per molecule, no atomics needed.
   * @param virial 9-double device global accumulator, atomicAdd'd
   *   (contended across molecules) -- always computed; whether it's
   *   *used* (added into conf.old().virial_tensor) is the host's
   *   decision, matching the CPU's `V == math::atomic_virial` gate.
   * @param error_flag device int, atomicExch'd to a nonzero code (1:
   *   orthogonal reference vectors, matching CPU's E_SHAKE_FAILURE; 2:
   *   exceeded max_iterations) -- never reset by this kernel, caller
   *   zeroes before launch and checks after.
   */
  void launch_shake_solvent(
      double3* pos,
      const double3* old_pos,
      const ShakeConstraint* constraints,
      unsigned num_constraints,
      const double* inv_mass_local,
      unsigned num_atoms_per_molecule,
      unsigned first_atom,
      unsigned num_molecules,
      double tolerance,
      unsigned max_iterations,
      math::boundary_enum boundary,
      math::Box box,
      double dt2,
      double3* constraint_force,
      double* virial,
      int* error_flag,
      cudaStream_t stream = 0);

} // namespace gpu
