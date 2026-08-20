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
 * @file position_restraint_kernels.h
 * Host-callable entry point for the position-restraint force/energy
 * kernel. One thread per restraint *term* (topo.position_restraints(),
 * a sparse list of restrained atoms -- not one thread per atom), direct
 * global atomicAdd -- same shape as angle_kernels.h/quartic_bond_kernels.h
 * (small, static term list; no pairlist).
 *
 * Exact CPU formula (position_restraint_interaction.cc's
 * _calculate_position_restraint_interactions): for restraint term idx
 * restraining atom seq[idx],
 * v = nearest_image(pos(seq), ref[idx]),
 * f = -(force_constant / inv_bf_scale[idx]) * v,
 * force(seq) += f,
 * e = 0.5 * force_constant / inv_bf_scale[idx] * |v|^2, accumulated into
 * posrest_energy[atom_energy_group[idx]]. `ref`/`inv_bf_scale` are
 * compact, one entry per restraint term (not per atom) -- the caller
 * (CUDA_Position_Restraint_Interaction::init()) already resolved
 * conf.special().reference_positions/bfactors down to this term's index,
 * same as m_atom_energy_group already does for every other bonded-term
 * kernel in this directory. `inv_bf_scale[idx]` is 1.0 for every term
 * when posrest.posrest != posrest_bfactor (uploaded that way once in
 * init(), not branched on per-kernel-call) -- see that class's doc
 * comment.
 *
 * No virial contribution at all -- the CPU code has it explicitly
 * commented out ("should there be a contribution of this special ia to
 * the virial? no there should be NO contribution"), so this kernel never
 * touches a virial buffer.
 */

#pragma once

#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "math/gmath.h"

namespace gpu {

  void launch_position_restraint(
      math::CuVArray::View pos,
      const unsigned* seq,
      const FPL3_TYPE* ref,
      const FPL_TYPE* inv_bf_scale,
      FPL_TYPE force_constant,
      const unsigned* atom_energy_group,
      unsigned num_restraints,
      math::boundary_enum boundary,
      math::Box box,
      FPL3_TYPE* force,
      double* posrest_energy,
      cudaStream_t stream = 0);

} // namespace gpu
