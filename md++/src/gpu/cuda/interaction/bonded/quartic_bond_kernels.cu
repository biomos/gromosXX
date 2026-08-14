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
 * @file quartic_bond_kernels.cu
 * Quartic bond force/energy kernel. See quartic_bond_kernels.h.
 */

#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/memory/types.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/math/periodicity.h"
#include "math/gmath.h"

#include "quartic_bond_kernels.h"

namespace gpu {

template <math::boundary_enum BOUNDARY>
__global__ void quartic_bond_kernel(
    math::CuVArray::View pos,
    const unsigned* __restrict__ bond_i,
    const unsigned* __restrict__ bond_j,
    const unsigned* __restrict__ bond_type,
    const FPL_TYPE* __restrict__ K,
    const FPL_TYPE* __restrict__ r0,
    const unsigned* __restrict__ atom_energy_group,
    unsigned num_bonds,
    gpu::Periodicity<BOUNDARY> periodicity,
    FPL3_TYPE* force,
    double* bond_energy,
    double* virial) {

  const unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= num_bonds) return;

  const unsigned i = bond_i[idx];
  const unsigned j = bond_j[idx];
  const unsigned type = bond_type[idx];

  const FPL3_TYPE v = periodicity.nearest_image(pos(i), pos(j));

  // dist2 - r02 is a near-cancellation for a bond close to its
  // equilibrium length (dist2 ~ r02) -- compute it in double regardless
  // of FPL_TYPE (float under the default FP_PRECISION=1 build) so the
  // subtraction itself doesn't amplify the position mirror's float
  // truncation into a visibly wrong force; found via quartic_bond_gpu.
  // t.cc's CPU-vs-GPU comparison exceeding tolerance with FPL_TYPE-only
  // arithmetic here.
  const double dist2 = static_cast<double>(v.x) * v.x +
                        static_cast<double>(v.y) * v.y +
                        static_cast<double>(v.z) * v.z;
  const double r0d = static_cast<double>(r0[type]);
  const double r02 = r0d * r0d;
  const double delta = dist2 - r02;

  const double coeff = -static_cast<double>(K[type]) * delta;
  const FPL3_TYPE f = static_cast<FPL_TYPE>(coeff) * v;

  atomicAdd(&force[i].x, f.x);
  atomicAdd(&force[i].y, f.y);
  atomicAdd(&force[i].z, f.z);
  atomicAdd(&force[j].x, -f.x);
  atomicAdd(&force[j].y, -f.y);
  atomicAdd(&force[j].z, -f.z);

  const double e = 0.25 * static_cast<double>(K[type]) * delta * delta;
  atomicAdd(&bond_energy[atom_energy_group[i]], e);

  // Unconditional atomic virial, exact CPU formula (quartic_bond_interaction.cc):
  // virial(a, c) += v(a) * f(c).
  atomicAdd(&virial[0], static_cast<double>(v.x * f.x)); // (0,0)
  atomicAdd(&virial[1], static_cast<double>(v.x * f.y)); // (0,1)
  atomicAdd(&virial[2], static_cast<double>(v.x * f.z)); // (0,2)
  atomicAdd(&virial[3], static_cast<double>(v.y * f.x)); // (1,0)
  atomicAdd(&virial[4], static_cast<double>(v.y * f.y)); // (1,1)
  atomicAdd(&virial[5], static_cast<double>(v.y * f.z)); // (1,2)
  atomicAdd(&virial[6], static_cast<double>(v.z * f.x)); // (2,0)
  atomicAdd(&virial[7], static_cast<double>(v.z * f.y)); // (2,1)
  atomicAdd(&virial[8], static_cast<double>(v.z * f.z)); // (2,2)
}

} // namespace gpu

void gpu::launch_quartic_bond(
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
    FPL3_TYPE* force,
    double* bond_energy,
    double* virial,
    cudaStream_t stream) {

  if (num_bonds == 0) return;

  const unsigned threads = 256;
  const unsigned blocks = (num_bonds + threads - 1) / threads;

  switch (boundary) {
    case math::vacuum:
      gpu::quartic_bond_kernel<math::vacuum><<<blocks, threads, 0, stream>>>(
          pos, bond_i, bond_j, bond_type, K, r0, atom_energy_group, num_bonds,
          gpu::Periodicity<math::vacuum>(box), force, bond_energy, virial);
      break;
    case math::rectangular:
      gpu::quartic_bond_kernel<math::rectangular><<<blocks, threads, 0, stream>>>(
          pos, bond_i, bond_j, bond_type, K, r0, atom_energy_group, num_bonds,
          gpu::Periodicity<math::rectangular>(box), force, bond_energy, virial);
      break;
    case math::triclinic:
      gpu::quartic_bond_kernel<math::triclinic><<<blocks, threads, 0, stream>>>(
          pos, bond_i, bond_j, bond_type, K, r0, atom_energy_group, num_bonds,
          gpu::Periodicity<math::triclinic>(box), force, bond_energy, virial);
      break;
    default:
      break;
  }
}
