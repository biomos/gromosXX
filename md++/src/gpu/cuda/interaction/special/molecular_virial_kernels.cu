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
 * @file molecular_virial_kernels.cu
 * Molecular-virial correction kernels. See molecular_virial_kernels.h.
 */

#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/memory/types.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/math/periodicity.h"
#include "math/gmath.h"

#include "molecular_virial_kernels.h"

namespace gpu {

template <math::boundary_enum BOUNDARY>
__global__ void group_com_kernel(
    math::CuVArray::View pos,
    const float* __restrict__ mass,
    const unsigned* __restrict__ group_offsets,
    unsigned num_groups,
    gpu::Periodicity<BOUNDARY> periodicity,
    FPL3_TYPE* com_pos) {

  const unsigned g = blockIdx.x * blockDim.x + threadIdx.x;
  if (g >= num_groups) return;

  const unsigned start = group_offsets[g];
  const unsigned end = group_offsets[g + 1];

  FPL3_TYPE prev = pos(start);
  double com_x = 0.0, com_y = 0.0, com_z = 0.0, tot_mass = 0.0;

  for (unsigned a = start; a < end; ++a) {
    const double m = static_cast<double>(mass[a]);
    tot_mass += m;

    const FPL3_TYPE p = periodicity.nearest_image(pos(a), prev);

    com_x += m * (static_cast<double>(p.x) + static_cast<double>(prev.x));
    com_y += m * (static_cast<double>(p.y) + static_cast<double>(prev.y));
    com_z += m * (static_cast<double>(p.z) + static_cast<double>(prev.z));

    prev.x += p.x;
    prev.y += p.y;
    prev.z += p.z;
  }

  com_pos[g].x = static_cast<FPL_TYPE>(com_x / tot_mass);
  com_pos[g].y = static_cast<FPL_TYPE>(com_y / tot_mass);
  com_pos[g].z = static_cast<FPL_TYPE>(com_z / tot_mass);
}

template <math::boundary_enum BOUNDARY>
__global__ void molecular_virial_correction_kernel(
    math::CuVArray::View pos,
    const FPH3_TYPE* __restrict__ force,
    const unsigned* __restrict__ group_id,
    const FPL3_TYPE* __restrict__ com_pos,
    unsigned num_atoms,
    gpu::Periodicity<BOUNDARY> periodicity,
    double* virial) {

  const unsigned a = blockIdx.x * blockDim.x + threadIdx.x;
  if (a >= num_atoms) return;

  const unsigned g = group_id[a];
  const FPL3_TYPE r = periodicity.nearest_image(pos(a), com_pos[g]);
  const FPH3_TYPE f = force[a];

  const double rx = r.x, ry = r.y, rz = r.z;
  const double fx = f.x, fy = f.y, fz = f.z;

  // Exact CPU formula (util/prepare_virial.cc's _atomic_to_molecular_
  // virial()): corrP(b, a_component) += force(a)(a_component) * r(b),
  // then virial_tensor -= corrP -- accumulate the already-negated value
  // directly (row-major b*3+a_component, matching virial_accumulate_
  // kernels.h's convention) so the caller can feed this straight into
  // launch_accumulate_virial9().
  atomicAdd(&virial[0], -(rx * fx)); // (0,0)
  atomicAdd(&virial[1], -(rx * fy)); // (0,1)
  atomicAdd(&virial[2], -(rx * fz)); // (0,2)
  atomicAdd(&virial[3], -(ry * fx)); // (1,0)
  atomicAdd(&virial[4], -(ry * fy)); // (1,1)
  atomicAdd(&virial[5], -(ry * fz)); // (1,2)
  atomicAdd(&virial[6], -(rz * fx)); // (2,0)
  atomicAdd(&virial[7], -(rz * fy)); // (2,1)
  atomicAdd(&virial[8], -(rz * fz)); // (2,2)
}

} // namespace gpu

void gpu::launch_group_com(
    math::CuVArray::View pos,
    const float* mass,
    const unsigned* group_offsets,
    unsigned num_groups,
    math::boundary_enum boundary,
    math::Box box,
    FPL3_TYPE* com_pos,
    cudaStream_t stream) {

  if (num_groups == 0) return;

  const unsigned threads = 128;
  const unsigned blocks = (num_groups + threads - 1) / threads;

  switch (boundary) {
    case math::vacuum:
      gpu::group_com_kernel<math::vacuum><<<blocks, threads, 0, stream>>>(
          pos, mass, group_offsets, num_groups, gpu::Periodicity<math::vacuum>(box), com_pos);
      break;
    case math::rectangular:
      gpu::group_com_kernel<math::rectangular><<<blocks, threads, 0, stream>>>(
          pos, mass, group_offsets, num_groups, gpu::Periodicity<math::rectangular>(box), com_pos);
      break;
    case math::triclinic:
      gpu::group_com_kernel<math::triclinic><<<blocks, threads, 0, stream>>>(
          pos, mass, group_offsets, num_groups, gpu::Periodicity<math::triclinic>(box), com_pos);
      break;
    default:
      break;
  }
}

void gpu::launch_molecular_virial_correction(
    math::CuVArray::View pos,
    const FPH3_TYPE* force,
    const unsigned* group_id,
    const FPL3_TYPE* com_pos,
    unsigned num_atoms,
    math::boundary_enum boundary,
    math::Box box,
    double* virial,
    cudaStream_t stream) {

  if (num_atoms == 0) return;

  const unsigned threads = 256;
  const unsigned blocks = (num_atoms + threads - 1) / threads;

  switch (boundary) {
    case math::vacuum:
      gpu::molecular_virial_correction_kernel<math::vacuum><<<blocks, threads, 0, stream>>>(
          pos, force, group_id, com_pos, num_atoms, gpu::Periodicity<math::vacuum>(box), virial);
      break;
    case math::rectangular:
      gpu::molecular_virial_correction_kernel<math::rectangular><<<blocks, threads, 0, stream>>>(
          pos, force, group_id, com_pos, num_atoms, gpu::Periodicity<math::rectangular>(box), virial);
      break;
    case math::triclinic:
      gpu::molecular_virial_correction_kernel<math::triclinic><<<blocks, threads, 0, stream>>>(
          pos, force, group_id, com_pos, num_atoms, gpu::Periodicity<math::triclinic>(box), virial);
      break;
    default:
      break;
  }
}
