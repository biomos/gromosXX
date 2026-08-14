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
 * @file angle_kernels.cu
 * Harmonic (cosine) bond-angle force/energy kernel. See angle_kernels.h.
 */

#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/memory/types.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/math/periodicity.h"
#include "math/gmath.h"

#include "angle_kernels.h"

namespace gpu {

template <math::boundary_enum BOUNDARY>
__global__ void angle_kernel(
    math::CuVArray::View pos,
    const unsigned* __restrict__ angle_i,
    const unsigned* __restrict__ angle_j,
    const unsigned* __restrict__ angle_k,
    const unsigned* __restrict__ angle_type,
    const FPL_TYPE* __restrict__ K,
    const FPL_TYPE* __restrict__ cos0,
    const unsigned* __restrict__ atom_energy_group,
    unsigned num_angles,
    gpu::Periodicity<BOUNDARY> periodicity,
    FPL3_TYPE* force,
    double* angle_energy,
    double* virial) {

  const unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= num_angles) return;

  const unsigned i = angle_i[idx];
  const unsigned j = angle_j[idx];
  const unsigned k = angle_k[idx];
  const unsigned type = angle_type[idx];

  const FPL3_TYPE rij = periodicity.nearest_image(pos(i), pos(j));
  const FPL3_TYPE rkj = periodicity.nearest_image(pos(k), pos(j));

  // Scalar geometry (dij, dkj, cost) done in double -- cost - cos0 is a
  // near-cancellation for an angle close to its equilibrium value, same
  // rationale as quartic_bond_kernels.cu's dist2 - r0^2 (see that file's
  // doc comment for the confirmed magnitude of the resulting error under
  // the default FP_PRECISION=1 float position mirror).
  const double rij_x = rij.x, rij_y = rij.y, rij_z = rij.z;
  const double rkj_x = rkj.x, rkj_y = rkj.y, rkj_z = rkj.z;

  const double dij2 = rij_x*rij_x + rij_y*rij_y + rij_z*rij_z;
  const double dkj2 = rkj_x*rkj_x + rkj_y*rkj_y + rkj_z*rkj_z;
  const double dij = sqrt(dij2);
  const double dkj = sqrt(dkj2);

  const double ip = rij_x*rkj_x + rij_y*rkj_y + rij_z*rkj_z;
  const double cost = ip / (dij * dkj);

  const double Kd = static_cast<double>(K[type]);
  const double cos0d = static_cast<double>(cos0[type]);
  const double df = -Kd * (cost - cos0d);

  const double fi_x = (df/dij) * (rkj_x/dkj - rij_x/dij*cost);
  const double fi_y = (df/dij) * (rkj_y/dkj - rij_y/dij*cost);
  const double fi_z = (df/dij) * (rkj_z/dkj - rij_z/dij*cost);

  const double fk_x = (df/dkj) * (rij_x/dij - rkj_x/dkj*cost);
  const double fk_y = (df/dkj) * (rij_y/dij - rkj_y/dkj*cost);
  const double fk_z = (df/dkj) * (rij_z/dij - rkj_z/dkj*cost);

  const double fj_x = -fi_x - fk_x;
  const double fj_y = -fi_y - fk_y;
  const double fj_z = -fi_z - fk_z;

  atomicAdd(&force[i].x, static_cast<FPL_TYPE>(fi_x));
  atomicAdd(&force[i].y, static_cast<FPL_TYPE>(fi_y));
  atomicAdd(&force[i].z, static_cast<FPL_TYPE>(fi_z));
  atomicAdd(&force[j].x, static_cast<FPL_TYPE>(fj_x));
  atomicAdd(&force[j].y, static_cast<FPL_TYPE>(fj_y));
  atomicAdd(&force[j].z, static_cast<FPL_TYPE>(fj_z));
  atomicAdd(&force[k].x, static_cast<FPL_TYPE>(fk_x));
  atomicAdd(&force[k].y, static_cast<FPL_TYPE>(fk_y));
  atomicAdd(&force[k].z, static_cast<FPL_TYPE>(fk_z));

  const double delta = cost - cos0d;
  const double e = 0.5 * Kd * delta * delta;
  atomicAdd(&angle_energy[atom_energy_group[idx]], e);

  // Unconditional atomic virial, exact CPU formula (angle_interaction.cc):
  // virial(a, b) += rij(a)*fi(b) + rkj(a)*fk(b).
  atomicAdd(&virial[0], rij_x*fi_x + rkj_x*fk_x); // (0,0)
  atomicAdd(&virial[1], rij_x*fi_y + rkj_x*fk_y); // (0,1)
  atomicAdd(&virial[2], rij_x*fi_z + rkj_x*fk_z); // (0,2)
  atomicAdd(&virial[3], rij_y*fi_x + rkj_y*fk_x); // (1,0)
  atomicAdd(&virial[4], rij_y*fi_y + rkj_y*fk_y); // (1,1)
  atomicAdd(&virial[5], rij_y*fi_z + rkj_y*fk_z); // (1,2)
  atomicAdd(&virial[6], rij_z*fi_x + rkj_z*fk_x); // (2,0)
  atomicAdd(&virial[7], rij_z*fi_y + rkj_z*fk_y); // (2,1)
  atomicAdd(&virial[8], rij_z*fi_z + rkj_z*fk_z); // (2,2)
}

} // namespace gpu

void gpu::launch_angle(
    math::CuVArray::View pos,
    const unsigned* angle_i,
    const unsigned* angle_j,
    const unsigned* angle_k,
    const unsigned* angle_type,
    const FPL_TYPE* K,
    const FPL_TYPE* cos0,
    const unsigned* atom_energy_group,
    unsigned num_angles,
    math::boundary_enum boundary,
    math::Box box,
    FPL3_TYPE* force,
    double* angle_energy,
    double* virial,
    cudaStream_t stream) {

  if (num_angles == 0) return;

  const unsigned threads = 256;
  const unsigned blocks = (num_angles + threads - 1) / threads;

  switch (boundary) {
    case math::vacuum:
      gpu::angle_kernel<math::vacuum><<<blocks, threads, 0, stream>>>(
          pos, angle_i, angle_j, angle_k, angle_type, K, cos0, atom_energy_group,
          num_angles, gpu::Periodicity<math::vacuum>(box), force, angle_energy, virial);
      break;
    case math::rectangular:
      gpu::angle_kernel<math::rectangular><<<blocks, threads, 0, stream>>>(
          pos, angle_i, angle_j, angle_k, angle_type, K, cos0, atom_energy_group,
          num_angles, gpu::Periodicity<math::rectangular>(box), force, angle_energy, virial);
      break;
    case math::triclinic:
      gpu::angle_kernel<math::triclinic><<<blocks, threads, 0, stream>>>(
          pos, angle_i, angle_j, angle_k, angle_type, K, cos0, atom_energy_group,
          num_angles, gpu::Periodicity<math::triclinic>(box), force, angle_energy, virial);
      break;
    default:
      break;
  }
}
