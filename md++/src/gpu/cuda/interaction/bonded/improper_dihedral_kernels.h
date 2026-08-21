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
 * @file improper_dihedral_kernels.h
 * Host-callable entry point for the improper dihedral force/energy
 * kernel (PLAN.md §10 step 14, bonded forces). One thread per term,
 * direct global atomicAdd -- same shape as quartic_bond_kernels.h/
 * angle_kernels.h.
 *
 * Exact CPU formula (improper_dihedral_interaction.cc's
 * _calculate_improper_interactions) -- see that file for the full
 * derivation; not restated here. Two behavioural differences from the
 * CPU path, both because a `__global__` kernel cannot call
 * `io::messages.add()`:
 *  - `acs` (the cosine of the improper angle) is clamped to [-1, 1] with
 *    no `io::message::critical` if it strays past `1 + math::epsilon` --
 *    the CPU path treats that as a fatal error; this kernel just clamps.
 *    Not expected to trigger for any topology this port is validated
 *    against, but worth knowing if a future topology hits it.
 *  - the "bond angle close to 180 degrees" `ki`/`kl` zeroing (dmj2/dnk2
 *    below `1e-10 * dkj2`) is applied silently, without the CPU path's
 *    `io::message::warning`.
 */

#pragma once

#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "math/gmath.h"

namespace gpu {

  void launch_improper_dihedral(
      math::CuVArray::View pos,
      const unsigned* dihedral_i,
      const unsigned* dihedral_j,
      const unsigned* dihedral_k,
      const unsigned* dihedral_l,
      const unsigned* dihedral_type,
      const FPL_TYPE* K,
      const FPL_TYPE* q0,
      const unsigned* atom_energy_group,
      unsigned num_dihedrals,
      math::boundary_enum boundary,
      math::Box box,
      FPH3_TYPE* force,
      double* improper_energy,
      double* virial,
      cudaStream_t stream = 0);

} // namespace gpu
