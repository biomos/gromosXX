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
 * @file improper_dihedral_kernels.cu
 * Improper dihedral force/energy kernel. See improper_dihedral_kernels.h.
 */

#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/memory/types.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/math/periodicity.h"
#include "math/gmath.h"

#include "improper_dihedral_kernels.h"

namespace gpu {

template <math::boundary_enum BOUNDARY>
__global__ void improper_dihedral_kernel(
    math::CuVArray::View pos,
    const unsigned* __restrict__ dihedral_i,
    const unsigned* __restrict__ dihedral_j,
    const unsigned* __restrict__ dihedral_k,
    const unsigned* __restrict__ dihedral_l,
    const unsigned* __restrict__ dihedral_type,
    const FPL_TYPE* __restrict__ K,
    const FPL_TYPE* __restrict__ q0,
    const unsigned* __restrict__ atom_energy_group,
    unsigned num_dihedrals,
    gpu::Periodicity<BOUNDARY> periodicity,
    FPH3_TYPE* force,
    double* improper_energy,
    double* virial) {

  const unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= num_dihedrals) return;

  const unsigned i = dihedral_i[idx];
  const unsigned j = dihedral_j[idx];
  const unsigned k = dihedral_k[idx];
  const unsigned l = dihedral_l[idx];
  const unsigned type = dihedral_type[idx];

  // All geometry in double -- q - q0 (via acos) is sensitive near the
  // equilibrium improper angle, same rationale as quartic_bond_kernels.
  // cu/angle_kernels.cu.
  const FPL3_TYPE rkj_f = periodicity.nearest_image(pos(k), pos(j));
  const FPL3_TYPE rij_f = periodicity.nearest_image(pos(i), pos(j));
  const FPL3_TYPE rkl_f = periodicity.nearest_image(pos(k), pos(l));
  const FPL3_TYPE rlj_f = periodicity.nearest_image(pos(l), pos(j));

  const double rkj_x = rkj_f.x, rkj_y = rkj_f.y, rkj_z = rkj_f.z;
  const double rij_x = rij_f.x, rij_y = rij_f.y, rij_z = rij_f.z;
  const double rkl_x = rkl_f.x, rkl_y = rkl_f.y, rkl_z = rkl_f.z;
  const double rlj_x = rlj_f.x, rlj_y = rlj_f.y, rlj_z = rlj_f.z;

  // rmj = cross(rij, rkj), rnk = cross(rkj, rkl)
  const double rmj_x = rij_y*rkj_z - rij_z*rkj_y;
  const double rmj_y = rij_z*rkj_x - rij_x*rkj_z;
  const double rmj_z = rij_x*rkj_y - rij_y*rkj_x;

  const double rnk_x = rkj_y*rkl_z - rkj_z*rkl_y;
  const double rnk_y = rkj_z*rkl_x - rkj_x*rkl_z;
  const double rnk_z = rkj_x*rkl_y - rkj_y*rkl_x;

  const double dkj2 = rkj_x*rkj_x + rkj_y*rkj_y + rkj_z*rkj_z;
  const double dmj2 = rmj_x*rmj_x + rmj_y*rmj_y + rmj_z*rmj_z;
  const double dnk2 = rnk_x*rnk_x + rnk_y*rnk_y + rnk_z*rnk_z;
  const double dkj = sqrt(dkj2);
  const double dmj = sqrt(dmj2);
  const double dnk = sqrt(dnk2);

  double ip = rmj_x*rnk_x + rmj_y*rnk_y + rmj_z*rnk_z;
  double acs = ip / (dmj * dnk);
  // Clamp, no io::messages (not callable from device) -- see this file's
  // header doc comment.
  if (acs > 1.0) acs = 1.0;
  if (acs < -1.0) acs = -1.0;

  double q = acos(acs);

  const double ip2 = rij_x*rnk_x + rij_y*rnk_y + rij_z*rnk_z;
  if (ip2 < 0) q *= -1.0;

  const double Kd = static_cast<double>(K[type]);
  const double q0d = static_cast<double>(q0[type]);
  const double delta = q - q0d;

  double ki = -Kd * delta * dkj;
  double kl = -ki;
  if (dmj2 < (1.0e-10 * dkj2)) {
    ki = 0.0;
  } else {
    ki = ki / dmj2;
  }
  if (dnk2 < (1.0e-10 * dkj2)) {
    kl = 0.0;
  } else {
    kl = kl / dnk2;
  }

  const double kj1 = (rij_x*rkj_x + rij_y*rkj_y + rij_z*rkj_z) / dkj2 - 1.0;
  const double kj2 = (rkl_x*rkj_x + rkl_y*rkj_y + rkl_z*rkj_z) / dkj2;

  const double fi_x = ki * rmj_x, fi_y = ki * rmj_y, fi_z = ki * rmj_z;
  const double fl_x = kl * rnk_x, fl_y = kl * rnk_y, fl_z = kl * rnk_z;
  const double fj_x = kj1*fi_x - kj2*fl_x;
  const double fj_y = kj1*fi_y - kj2*fl_y;
  const double fj_z = kj1*fi_z - kj2*fl_z;
  const double fk_x = -(fi_x + fj_x + fl_x);
  const double fk_y = -(fi_y + fj_y + fl_y);
  const double fk_z = -(fi_z + fj_z + fl_z);

  atomicAdd(&force[i].x, static_cast<FPH_TYPE>(fi_x));
  atomicAdd(&force[i].y, static_cast<FPH_TYPE>(fi_y));
  atomicAdd(&force[i].z, static_cast<FPH_TYPE>(fi_z));
  atomicAdd(&force[j].x, static_cast<FPH_TYPE>(fj_x));
  atomicAdd(&force[j].y, static_cast<FPH_TYPE>(fj_y));
  atomicAdd(&force[j].z, static_cast<FPH_TYPE>(fj_z));
  atomicAdd(&force[k].x, static_cast<FPH_TYPE>(fk_x));
  atomicAdd(&force[k].y, static_cast<FPH_TYPE>(fk_y));
  atomicAdd(&force[k].z, static_cast<FPH_TYPE>(fk_z));
  atomicAdd(&force[l].x, static_cast<FPH_TYPE>(fl_x));
  atomicAdd(&force[l].y, static_cast<FPH_TYPE>(fl_y));
  atomicAdd(&force[l].z, static_cast<FPH_TYPE>(fl_z));

  const double e = 0.5 * Kd * delta * delta;
  atomicAdd(&improper_energy[atom_energy_group[idx]], e);

  // Unconditional atomic virial, exact CPU formula (improper_dihedral_
  // interaction.cc): virial(a,b) += rij(a)*fi(b) + rkj(a)*fk(b) + rlj(a)*fl(b).
  atomicAdd(&virial[0], rij_x*fi_x + rkj_x*fk_x + rlj_x*fl_x); // (0,0)
  atomicAdd(&virial[1], rij_x*fi_y + rkj_x*fk_y + rlj_x*fl_y); // (0,1)
  atomicAdd(&virial[2], rij_x*fi_z + rkj_x*fk_z + rlj_x*fl_z); // (0,2)
  atomicAdd(&virial[3], rij_y*fi_x + rkj_y*fk_x + rlj_y*fl_x); // (1,0)
  atomicAdd(&virial[4], rij_y*fi_y + rkj_y*fk_y + rlj_y*fl_y); // (1,1)
  atomicAdd(&virial[5], rij_y*fi_z + rkj_y*fk_z + rlj_y*fl_z); // (1,2)
  atomicAdd(&virial[6], rij_z*fi_x + rkj_z*fk_x + rlj_z*fl_x); // (2,0)
  atomicAdd(&virial[7], rij_z*fi_y + rkj_z*fk_y + rlj_z*fl_y); // (2,1)
  atomicAdd(&virial[8], rij_z*fi_z + rkj_z*fk_z + rlj_z*fl_z); // (2,2)
}

} // namespace gpu

void gpu::launch_improper_dihedral(
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
    cudaStream_t stream) {

  if (num_dihedrals == 0) return;

  const unsigned threads = 256;
  const unsigned blocks = (num_dihedrals + threads - 1) / threads;

  switch (boundary) {
    case math::vacuum:
      gpu::improper_dihedral_kernel<math::vacuum><<<blocks, threads, 0, stream>>>(
          pos, dihedral_i, dihedral_j, dihedral_k, dihedral_l, dihedral_type, K, q0,
          atom_energy_group, num_dihedrals, gpu::Periodicity<math::vacuum>(box),
          force, improper_energy, virial);
      break;
    case math::rectangular:
      gpu::improper_dihedral_kernel<math::rectangular><<<blocks, threads, 0, stream>>>(
          pos, dihedral_i, dihedral_j, dihedral_k, dihedral_l, dihedral_type, K, q0,
          atom_energy_group, num_dihedrals, gpu::Periodicity<math::rectangular>(box),
          force, improper_energy, virial);
      break;
    case math::triclinic:
      gpu::improper_dihedral_kernel<math::triclinic><<<blocks, threads, 0, stream>>>(
          pos, dihedral_i, dihedral_j, dihedral_k, dihedral_l, dihedral_type, K, q0,
          atom_energy_group, num_dihedrals, gpu::Periodicity<math::triclinic>(box),
          force, improper_energy, virial);
      break;
    default:
      break;
  }
}
