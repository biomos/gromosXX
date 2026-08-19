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
 * @file m_shake_kernels.h
 * Host-callable entry point for the M-SHAKE solvent kernel -- matches
 * algorithm::M_Shake::m_shake_molecule()/solvent() (m_shake.cc)
 * exactly: solves all 3 constraints of a molecule simultaneously each
 * iteration via a direct 3x3 matrix inversion, instead of SHAKE's
 * one-constraint-at-a-time Gauss-Seidel sweep. Still iterative (the
 * A matrix depends on the current, not-yet-converged positions), but
 * converges in far fewer iterations than plain SHAKE since the whole
 * coupled system is solved each pass -- for exactly this reason (a
 * shorter per-thread dependency chain), this was ported as the
 * follow-up to CUDA_Shake's solvent path once profiling showed that
 * path's cost dominated by iteration count under a grid too small to
 * hide per-iteration global-memory latency (KNOWN_ISSUES.md).
 *
 * The `factor` 3x3 matrix, `mass_i`, and `constr_length2` are the same
 * for every molecule of the (single, GROMOS-required-identical)
 * solvent type -- computed once host-side in CUDA_M_Shake::init(),
 * uploaded once, and reused by every thread every step, exactly
 * mirroring the CPU class's own member variables.
 *
 * v1 scope matches algorithm::M_Shake::init()'s own gates (checked
 * host-side in CUDA_M_Shake::init(), not repeated here): exactly one
 * solvent type, exactly 3 atoms per molecule, exactly 3 distance
 * constraints.
 */

#pragma once

#include "gpu/cuda/memory/cuvector.h"
#include "math/gmath.h"

namespace gpu {

  struct MShakeConstraint {
    unsigned i;
    unsigned j;
  };

  /**
   * @param pos current positions, flattened double3, device, updated
   *   in place for the solvent range only.
   * @param old_pos reference (previous-step) positions, same layout,
   *   read-only.
   * @param constr the molecule-local (0/1/2) atom index pairs of the 3
   *   distance constraints, same order as `factor`/`constr_length2`'s
   *   rows/entries.
   * @param factor row-major 3x3, `factor(k,l)` at `factor[3*k+l]`,
   *   matches algorithm::M_Shake::factor exactly (see init()).
   * @param constr_length2 reference (squared) constraint lengths,
   *   indexed the same way as `constr`.
   * @param mass_i inverse masses of the molecule's 3 atoms, indexed by
   *   molecule-local atom index (0/1/2), not by constraint.
   * @param first_atom global index of molecule 0's first atom
   *   (== topo.num_solute_atoms()).
   * @param num_molecules number of (single) solvent-type molecules.
   * @param tolerance, max_iterations: algorithm::M_Shake's own
   *   `tolerance()`/`max_iterations()`.
   * @param dt2i 1 / (dt*dt).
   * @param do_virial matches the CPU's `pcouple.virial ==
   *   atomic_virial` gate, evaluated host-side.
   * @param constraint_force written into (device, global, solvent-
   *   range-sized; disjoint per molecule, no atomics needed) -- NOT
   *   pre-zeroed by this kernel, matching the CPU's own `+=`
   *   accumulation into a caller-zeroed buffer.
   * @param virial 9-double device global accumulator. Reduced once per
   *   block (shared memory) before a single atomicAdd per block per
   *   component, not once per molecule -- the naive per-molecule
   *   atomicAdd pattern this avoids was already found to be a real
   *   contention bottleneck in the sibling SHAKE/SETTLE kernels.
   * @param error_flag device int, atomicExch'd to 1 if any molecule's
   *   matrix is singular (matches the CPU's `E_SHAKE_FAILURE` check).
   */
  void launch_m_shake_solvent(
      double3* pos,
      const double3* old_pos,
      const gpu::MShakeConstraint* constr,
      const double* factor,
      const double* constr_length2,
      const double* mass_i,
      unsigned first_atom,
      unsigned num_molecules,
      double tolerance,
      unsigned max_iterations,
      double dt2i,
      bool do_virial,
      double3* constraint_force,
      double* virial,
      int* error_flag,
      cudaStream_t stream = 0);

} // namespace gpu
