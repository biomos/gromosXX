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
 * @file position_restraint_kernels.cu
 * Position-restraint force/energy kernel. See position_restraint_kernels.h.
 */

#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/memory/types.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/math/periodicity.h"
#include "math/gmath.h"

#include "position_restraint_kernels.h"

namespace gpu {

template <math::boundary_enum BOUNDARY>
__global__ void position_restraint_kernel(
    math::CuVArray::View pos,
    const unsigned* __restrict__ seq,
    const FPL3_TYPE* __restrict__ ref,
    const FPL_TYPE* __restrict__ inv_bf_scale,
    FPL_TYPE force_constant,
    const unsigned* __restrict__ atom_energy_group,
    unsigned num_restraints,
    gpu::Periodicity<BOUNDARY> periodicity,
    FPL3_TYPE* force,
    double* posrest_energy) {

  const unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= num_restraints) return;

  const unsigned a = seq[idx];

  const FPL3_TYPE v = periodicity.nearest_image(pos(a), ref[idx]);

  const double v_x = v.x, v_y = v.y, v_z = v.z;
  const double k = static_cast<double>(force_constant) * static_cast<double>(inv_bf_scale[idx]);

  const double f_x = -k * v_x;
  const double f_y = -k * v_y;
  const double f_z = -k * v_z;

  atomicAdd(&force[a].x, static_cast<FPL_TYPE>(f_x));
  atomicAdd(&force[a].y, static_cast<FPL_TYPE>(f_y));
  atomicAdd(&force[a].z, static_cast<FPL_TYPE>(f_z));

  const double e = 0.5 * k * (v_x*v_x + v_y*v_y + v_z*v_z);
  atomicAdd(&posrest_energy[atom_energy_group[idx]], e);

  // No virial contribution -- exact CPU convention, see this file's doc
  // comment.
}

} // namespace gpu

void gpu::launch_position_restraint(
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
    cudaStream_t stream) {

  if (num_restraints == 0) return;

  const unsigned threads = 256;
  const unsigned blocks = (num_restraints + threads - 1) / threads;

  switch (boundary) {
    case math::vacuum:
      gpu::position_restraint_kernel<math::vacuum><<<blocks, threads, 0, stream>>>(
          pos, seq, ref, inv_bf_scale, force_constant, atom_energy_group,
          num_restraints, gpu::Periodicity<math::vacuum>(box), force, posrest_energy);
      break;
    case math::rectangular:
      gpu::position_restraint_kernel<math::rectangular><<<blocks, threads, 0, stream>>>(
          pos, seq, ref, inv_bf_scale, force_constant, atom_energy_group,
          num_restraints, gpu::Periodicity<math::rectangular>(box), force, posrest_energy);
      break;
    case math::triclinic:
      gpu::position_restraint_kernel<math::triclinic><<<blocks, threads, 0, stream>>>(
          pos, seq, ref, inv_bf_scale, force_constant, atom_energy_group,
          num_restraints, gpu::Periodicity<math::triclinic>(box), force, posrest_energy);
      break;
    default:
      break;
  }
}
