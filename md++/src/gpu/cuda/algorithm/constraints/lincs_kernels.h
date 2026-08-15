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
 * @file lincs_kernels.h
 * Host-callable entry points for the LINCS kernels (PLAN.md §10 step
 * 19, constraints), matching `_lincs<B>`/`_solve_lincs`
 * (algorithm/constraints/lincs.cc) exactly.
 *
 * Unlike SHAKE, LINCS's own CPU recursion (`_solve_lincs`'s `rec` loop)
 * is *already* Jacobi-shaped, not Gauss-Seidel: each round only reads
 * the previous round's solution vector (`rhs[1-w]`) and writes a new
 * one (`rhs[w]`) for every constraint, with no within-round ordering
 * dependency. That means, unlike solute SHAKE, this port doesn't need
 * to *reformulate* the algorithm to parallelize it -- it's the same
 * iteration structure as the CPU, just executed by one GPU thread per
 * constraint per round instead of a CPU loop, and *is* bit-comparable
 * to the CPU per constraint (up to floating-point non-associativity).
 *
 * A single "group" covers either the whole solute system (one
 * instance, `num_instances = 1`) or one solvent type (`num_instances =
 * num_molecules of that type`, all sharing the same per-molecule
 * coupling structure `topo.solvent(s).lincs()`, since every molecule
 * of a type has identical topology) -- see cuda_lincs.h for how the
 * two are set up. `A[i].a[n]` (the CPU's per-constraint-pair
 * coefficient cache) isn't stored on the GPU: `coef` (mass-derived,
 * static) is uploaded once per group, but `dot(B(i), B(coupled))`
 * (position-dependent, changes every call) is recomputed on the fly
 * inside `launch_lincs_round` instead of cached in a separate kernel --
 * coupling degree is small, so this trades a few redundant dot
 * products for a much smaller amount of GPU memory traffic overall.
 *
 * Deliberately uses plain `double`/`double3` throughout, not
 * `FPL_TYPE`/`FPL3_TYPE` -- same rationale as shake_kernels.h/
 * settle_kernels.h.
 */

#pragma once

#include "gpu/cuda/memory/cuvector.h"
#include "math/gmath.h"

namespace gpu {

  /**
   * One entry per LINCS constraint within a group, atom indices local
   * to the group's instance (0-indexed within the solute system, or
   * within one solvent molecule).
   */
  struct LincsConstraint {
    unsigned i;
    unsigned j;
    double r0;
    double sdiag;
    double mass_i;
    double mass_j;
  };

  /**
   * Computes `B(i) = normalized(nearest_image(old_pos(i), old_pos(j)))`
   * for every constraint in every instance of the group -- one thread
   * per (instance, local constraint). `B` is sized
   * `num_instances * num_constr_per_instance`.
   */
  void launch_lincs_compute_b(
      const double3* old_pos,
      const gpu::LincsConstraint* constraints,
      unsigned num_constr_per_instance,
      unsigned num_instances,
      unsigned first_atom,
      unsigned atom_stride_per_instance,
      math::boundary_enum boundary,
      math::Box box,
      double3* B,
      cudaStream_t stream = 0);

  /**
   * Computes the initial right-hand side/solution (`rhs[0]`/`sol`) from
   * the current (not-yet-corrected) positions -- exactly `_lincs<B>`'s
   * first `rhs[0](i) = sdiag[i] * (dot(B(i), r) - r0)` loop. One thread
   * per (instance, local constraint).
   */
  void launch_lincs_init_rhs(
      const double3* pos,
      const gpu::LincsConstraint* constraints,
      unsigned num_constr_per_instance,
      unsigned num_instances,
      unsigned first_atom,
      unsigned atom_stride_per_instance,
      math::boundary_enum boundary,
      math::Box box,
      const double3* B,
      double* rhs,
      double* sol,
      cudaStream_t stream = 0);

  /**
   * Computes the rotational-lengthening-corrected right-hand side --
   * exactly `_lincs<B>`'s second `rhs[0](i) = sdiag[i] * (r0 - p)` loop
   * (`p = sqrt(2*r0^2 - |r|^2)`, `p = 0` -- and a diagnostic counter
   * incremented, matching the CPU's "too much rotation" message -- if
   * the value under the square root is negative). One thread per
   * (instance, local constraint).
   */
  void launch_lincs_rotation_rhs(
      const double3* pos,
      const gpu::LincsConstraint* constraints,
      unsigned num_constr_per_instance,
      unsigned num_instances,
      unsigned first_atom,
      unsigned atom_stride_per_instance,
      math::boundary_enum boundary,
      math::Box box,
      double* rhs,
      double* sol,
      int* rotation_count,
      cudaStream_t stream = 0);

  /**
   * One round of `_solve_lincs`'s `rec` loop: `rhs_out(i) = sum_n
   * coef[i][n] * dot(B(i), B(coupled_n)) * rhs_in(coupled_n)`,
   * `sol(i) += rhs_out(i)`. Ping-pong `rhs_in`/`rhs_out` across
   * `lincs_order` calls (caller swaps pointers between calls -- this
   * kernel never reads and writes the same buffer). `coupled_offset`/
   * `coupled_index`/`coupled_coef` are the per-type (not per-instance)
   * flattened CSR coupling lists, `coupled_offset` sized
   * `num_constr_per_instance + 1`.
   */
  void launch_lincs_round(
      const double3* B,
      const unsigned* coupled_offset,
      const unsigned* coupled_index,
      const double* coupled_coef,
      unsigned num_constr_per_instance,
      unsigned num_instances,
      const double* rhs_in,
      double* rhs_out,
      double* sol,
      cudaStream_t stream = 0);

  /**
   * Applies the accumulated `sol` to `pos` -- exactly `_solve_lincs`'s
   * position-update loop (`pos(i) -= B(i)*sdiag[i]/mass_i*sol(i)`,
   * `pos(j) += ...`). `atomicAdd`s into `pos` since constraints sharing
   * an atom (a constrained chain) write the same position concurrently.
   */
  void launch_lincs_apply(
      double3* pos,
      const gpu::LincsConstraint* constraints,
      const double3* B,
      const double* sol,
      unsigned num_constr_per_instance,
      unsigned num_instances,
      unsigned first_atom,
      unsigned atom_stride_per_instance,
      cudaStream_t stream = 0);

} // namespace gpu
