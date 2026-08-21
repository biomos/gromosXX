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
 * @file quartic_bond_kernels.h
 * Host-callable entry point for the quartic bond force/energy kernel
 * (PLAN.md §10 step 14, bonded forces). One thread per bond term, direct
 * global atomicAdd (no shared-memory bucketing) -- term counts are
 * modest (aladip: ~80 bonds), unlike the nonbonded tile kernel's
 * per-block-of-many-pairs shape.
 *
 * Same host/device split convention as lj_crf_tiles.h: this header does
 * not include gpu/cuda/math/periodicity.h or declare the __global__
 * kernel itself, so it stays safe to include from a plain (non-nvcc)
 * .cc translation unit (e.g. cuda_quartic_bond_interaction.cc, which
 * compiles under groalgorithm/grointeraction -- no CUDA language
 * enabled there).
 *
 * Exact CPU formula (quartic_bond_interaction.cc's
 * _calculate_quartic_bond_interactions): for bond term (i, j, type),
 * v = nearest_image(pos(i), pos(j)), dist2 = |v|^2, r02 = r0[type]^2,
 * f = v * (-K[type] * (dist2 - r02)), force(i) += f, force(j) -= f,
 * e = 0.25 * K[type] * (dist2 - r02)^2, accumulated into
 * bond_energy[atom_energy_group[i]]. Virial is unconditional (the CPU
 * code's `if (V == math::atomic_virial)` gate is dead/commented-out):
 * virial(a, c) += v(a) * f(c).
 */

#pragma once

#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "math/gmath.h"

namespace gpu {

  void launch_quartic_bond(
      math::CuVArray::View pos,
      const unsigned* bond_i,
      const unsigned* bond_j,
      const unsigned* bond_type,
      const FPL_TYPE* K,
      const FPL_TYPE* r0,
      const unsigned* atom_energy_group,
      unsigned num_bonds,
      math::boundary_enum boundary,
      math::Box box,
      FPH3_TYPE* force,
      double* bond_energy,
      double* virial,
      cudaStream_t stream = 0);

} // namespace gpu
