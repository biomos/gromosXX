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
 * Host-callable entry points for the SHAKE kernels (PLAN.md §10 step
 * 15/17, constraints).
 *
 * `launch_shake_solvent`: one CUDA thread per solvent molecule, doing
 * the *exact same* sequential Gauss-Seidel iteration (with the
 * skip_now/skip_next per-atom convergence-tracking optimization) as
 * the CPU reference (algorithm::Shake::shake_iteration /
 * algorithm::Shake::solvent, shake.h) -- molecules are independent of
 * each other (no shared atoms), so parallelizing *across* molecules
 * while keeping each molecule's own inner loop sequential (exactly like
 * the CPU) is both correct and embarrassingly parallel.
 *
 * `launch_shake_solute_round`/`launch_shake_solute_apply`: solute
 * SHAKE's CPU reference (algorithm::Shake::solute()) is a fundamentally
 * different problem -- one large in-place Gauss-Seidel sweep over
 * *every* solute distance constraint at once, with genuine cross-
 * constraint data dependencies within a single iteration (updating
 * atom i's position mid-sweep changes what the very next constraint in
 * the list sees). These two kernels instead implement a **Jacobi-style**
 * parallel constraint solve: one thread per constraint reads a fixed
 * snapshot of positions (not updated mid-round) and atomicAdd's its
 * correction into a per-atom delta buffer (`launch_shake_solute_round`);
 * a second kernel then applies every atom's accumulated delta to the
 * position snapshot in one pass (`launch_shake_solute_apply`), so
 * atoms shared by multiple constraints (e.g. a constrained chain)
 * still get every constraint's contribution, just all computed
 * against the same starting point rather than sequentially. This
 * converges to the same constrained manifold as Gauss-Seidel (both are
 * standard iterative constraint solvers) but is **not** bit-comparable
 * to the CPU's specific iteration path -- see cuda_shake.h's doc
 * comment.
 *
 * Uses `FPL_TYPE`/`FPL3_TYPE` (gpu/cuda/memory/precision.h) for
 * per-atom position/force data, matching the rest of the GPU pipeline
 * (e.g. lj_crf_tiles.cu) -- this was previously deliberately plain
 * `double`/`double3` (position corrections directly become the next
 * MD step's integration input, so precision loss here is more
 * consequential than in an energy comparison), but profiling this
 * consumer GPU (RTX 5060 Ti) found double-precision throughput
 * catastrophically throttled relative to a CPU core (~19-135x slower
 * per constraint call, versus only ~1.5x for the already-mixed-
 * precision nonbonded kernels) -- consumer/GeForce parts run FP64 at
 * a small fraction of FP32 rate, unlike datacenter GPUs. Switched to
 * `FPL_TYPE` deliberately, accepting the precision tradeoff, per this
 * session's explicit direction. The global virial accumulator stays
 * plain `double` (matching lj_crf_tiles.cu's own convention: bulk
 * per-thread math in low precision, the one cross-molecule reduction
 * target in high precision).
 */

#pragma once

#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/memory/precision.h"
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
    FPL_TYPE r0sq;
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
      FPL3_TYPE* pos,
      const FPL3_TYPE* old_pos,
      const ShakeConstraint* constraints,
      unsigned num_constraints,
      const FPL_TYPE* inv_mass_local,
      unsigned num_atoms_per_molecule,
      unsigned first_atom,
      unsigned num_molecules,
      FPL_TYPE tolerance,
      unsigned max_iterations,
      math::boundary_enum boundary,
      math::Box box,
      FPL_TYPE dt2,
      FPH3_TYPE* constraint_force,
      double* virial,
      int* error_flag,
      cudaStream_t stream = 0);

  /**
   * One round of the Jacobi-style solute constraint solve: one thread
   * per constraint, reading `pos` (the fixed snapshot for this round --
   * NOT concurrently modified by this kernel) and atomicAdd-ing its
   * correction, scaled by inverse mass, into `delta` (device, global,
   * num_solute_atoms-sized, NOT zeroed by this kernel -- caller zeroes
   * before launch). Also atomicAdd's into `constraint_force`/`virial`
   * (raw, undivided by dt2 -- same convention as launch_shake_solvent)
   * and atomicExch's `changed_flag` to 1 if any constraint exceeded
   * tolerance (mirrors the CPU's `convergence = false`).
   * `launch_shake_solute_apply` (below) must be called afterwards, on
   * the same stream, to actually add `delta` into `pos` -- kept as two
   * kernels specifically so every constraint in this round sees the
   * same starting positions, matching the Jacobi (not Gauss-Seidel)
   * scheme this class deliberately uses (see this file's doc comment).
   */
  void launch_shake_solute_round(
      const FPL3_TYPE* pos,
      const FPL3_TYPE* old_pos,
      const ShakeConstraint* constraints,
      unsigned num_constraints,
      const FPL_TYPE* inv_mass,
      FPL_TYPE tolerance,
      math::boundary_enum boundary,
      math::Box box,
      FPL_TYPE dt2,
      FPL3_TYPE* delta,
      FPH3_TYPE* constraint_force,
      double* virial,
      int* changed_flag,
      int* error_flag,
      cudaStream_t stream = 0);

  /**
   * Applies `delta` (accumulated by launch_shake_solute_round) into
   * `pos` and resets `delta` to zero, one thread per solute atom.
   */
  void launch_shake_solute_apply(
      FPL3_TYPE* pos,
      FPL3_TYPE* delta,
      unsigned num_atoms,
      cudaStream_t stream = 0);

} // namespace gpu
