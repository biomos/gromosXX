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
 * @file dihedral_kernels.h
 * Host-callable entry point for the torsional dihedral (any m, any
 * shift angle) force/energy kernel (PLAN.md §10 step 14, bonded
 * forces), matching Dihedral_new_Interaction
 * (dihedral_new_interaction.cc). One thread per term, direct global
 * atomicAdd -- same shape as the other bonded-term kernels in this
 * directory.
 *
 * Two behavioural differences from the CPU path, both because a
 * `__global__` kernel cannot call `io::messages.add()`:
 *  - the "bond angle close to 180 degrees" `fi`/`fl` zeroing (dmj2/dnk2
 *    below `1e-10 * dkj2`) is applied silently, without the CPU path's
 *    `io::message::warning`.
 *  - dihedral-angle-minimum monitoring (`sim.param().print.
 *    monitor_dihedrals`, `conf.special().dihangle_trans`) is not
 *    implemented at all -- out of scope, hard-errored in
 *    CUDA_Dihedral_Interaction::init() if requested, same convention as
 *    every other CUDA gate in this codebase.
 */

#pragma once

#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "math/gmath.h"

namespace gpu {

  void launch_dihedral(
      math::CuVArray::View pos,
      const unsigned* dihedral_i,
      const unsigned* dihedral_j,
      const unsigned* dihedral_k,
      const unsigned* dihedral_l,
      const unsigned* dihedral_type,
      const FPL_TYPE* K,
      const FPL_TYPE* cospd,
      const FPL_TYPE* pd,
      const int* m,
      const unsigned* atom_energy_group,
      unsigned num_dihedrals,
      math::boundary_enum boundary,
      math::Box box,
      FPL3_TYPE* force,
      double* dihedral_energy,
      double* virial,
      cudaStream_t stream = 0);

} // namespace gpu
