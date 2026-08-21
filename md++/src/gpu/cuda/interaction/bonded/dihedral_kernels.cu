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
 * @file dihedral_kernels.cu
 * Torsional dihedral force/energy kernel. See dihedral_kernels.h.
 */

#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/memory/types.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/math/periodicity.h"
#include "math/gmath.h"

#include "dihedral_kernels.h"

namespace gpu {

template <math::boundary_enum BOUNDARY>
__global__ void dihedral_kernel(
    math::CuVArray::View pos,
    const unsigned* __restrict__ dihedral_i,
    const unsigned* __restrict__ dihedral_j,
    const unsigned* __restrict__ dihedral_k,
    const unsigned* __restrict__ dihedral_l,
    const unsigned* __restrict__ dihedral_type,
    const FPL_TYPE* __restrict__ K,
    const FPL_TYPE* __restrict__ cospd,
    const FPL_TYPE* __restrict__ pd,
    const int* __restrict__ m,
    const unsigned* __restrict__ atom_energy_group,
    unsigned num_dihedrals,
    gpu::Periodicity<BOUNDARY> periodicity,
    FPH3_TYPE* force,
    double* dihedral_energy,
    double* virial) {

  const unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= num_dihedrals) return;

  const unsigned i = dihedral_i[idx];
  const unsigned j = dihedral_j[idx];
  const unsigned k = dihedral_k[idx];
  const unsigned l = dihedral_l[idx];
  const unsigned type = dihedral_type[idx];

  const FPL3_TYPE rij_f = periodicity.nearest_image(pos(i), pos(j));
  const FPL3_TYPE rkj_f = periodicity.nearest_image(pos(k), pos(j));
  const FPL3_TYPE rkl_f = periodicity.nearest_image(pos(k), pos(l));
  const FPL3_TYPE rlj_f = periodicity.nearest_image(pos(l), pos(j));

  const double rij_x = rij_f.x, rij_y = rij_f.y, rij_z = rij_f.z;
  const double rkj_x = rkj_f.x, rkj_y = rkj_f.y, rkj_z = rkj_f.z;
  const double rkl_x = rkl_f.x, rkl_y = rkl_f.y, rkl_z = rkl_f.z;
  const double rlj_x = rlj_f.x, rlj_y = rlj_f.y, rlj_z = rlj_f.z;

  // rmj = cross(rij, rkj), rnk = cross(rkj, rkl)
  const double rmj_x = rij_y*rkj_z - rij_z*rkj_y;
  const double rmj_y = rij_z*rkj_x - rij_x*rkj_z;
  const double rmj_z = rij_x*rkj_y - rij_y*rkj_x;

  const double rnk_x = rkj_y*rkl_z - rkj_z*rkl_y;
  const double rnk_y = rkj_z*rkl_x - rkj_x*rkl_z;
  const double rnk_z = rkj_x*rkl_y - rkj_y*rkl_x;

  const double dmj2 = rmj_x*rmj_x + rmj_y*rmj_y + rmj_z*rmj_z;
  const double dnk2 = rnk_x*rnk_x + rnk_y*rnk_y + rnk_z*rnk_z;
  const double dkj2 = rkj_x*rkj_x + rkj_y*rkj_y + rkj_z*rkj_z;
  const double dkj = sqrt(dkj2);

  const double frim = (rij_x*rkj_x + rij_y*rkj_y + rij_z*rkj_z) / dkj2;
  const double frln = (rkl_x*rkj_x + rkl_y*rkj_y + rkl_z*rkj_z) / dkj2;

  const double rim_x = rij_x - frim*rkj_x;
  const double rim_y = rij_y - frim*rkj_y;
  const double rim_z = rij_z - frim*rkj_z;

  const double rln_x = frln*rkj_x - rkl_x;
  const double rln_y = frln*rkj_y - rkl_y;
  const double rln_z = frln*rkj_z - rkl_z;

  const double dim = sqrt(rim_x*rim_x + rim_y*rim_y + rim_z*rim_z);
  const double dln = sqrt(rln_x*rln_x + rln_y*rln_y + rln_z*rln_z);

  const double ip = rim_x*rln_x + rim_y*rln_y + rim_z*rln_z;
  double cosphi = ip / (dim * dln);
  if (cosphi > 1.0) cosphi = 1.0;
  if (cosphi < -1.0) cosphi = -1.0;
  double phi = acos(cosphi);

  const double sign = rij_x*rnk_x + rij_y*rnk_y + rij_z*rnk_z;
  if (sign < 0) phi *= -1.0;

  const double Kd = static_cast<double>(K[type]);
  const double deltad = static_cast<double>(pd[type]);
  const double md = static_cast<double>(m[type]);

  const double kj1 = frim - 1.0;
  const double kj2 = frln;
  const double ki = Kd * md * sin(md * phi - deltad);
  const double kl = -ki;

  double fi_x = ki * dkj * rmj_x;
  double fi_y = ki * dkj * rmj_y;
  double fi_z = ki * dkj * rmj_z;
  if (dmj2 < (1.0e-10 * dkj2)) {
    fi_x = fi_y = fi_z = 0.0;
  } else {
    fi_x /= dmj2; fi_y /= dmj2; fi_z /= dmj2;
  }

  double fl_x = kl * dkj * rnk_x;
  double fl_y = kl * dkj * rnk_y;
  double fl_z = kl * dkj * rnk_z;
  if (dnk2 < (1.0e-10 * dkj2)) {
    fl_x = fl_y = fl_z = 0.0;
  } else {
    fl_x /= dnk2; fl_y /= dnk2; fl_z /= dnk2;
  }

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

  const double e = Kd * (1.0 + cos(md * phi - deltad));
  atomicAdd(&dihedral_energy[atom_energy_group[idx]], e);

  // Unconditional atomic virial, exact CPU formula (dihedral_new_
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

void gpu::launch_dihedral(
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
    FPH3_TYPE* force,
    double* dihedral_energy,
    double* virial,
    cudaStream_t stream) {

  if (num_dihedrals == 0) return;

  const unsigned threads = 256;
  const unsigned blocks = (num_dihedrals + threads - 1) / threads;

  switch (boundary) {
    case math::vacuum:
      gpu::dihedral_kernel<math::vacuum><<<blocks, threads, 0, stream>>>(
          pos, dihedral_i, dihedral_j, dihedral_k, dihedral_l, dihedral_type,
          K, cospd, pd, m, atom_energy_group, num_dihedrals,
          gpu::Periodicity<math::vacuum>(box), force, dihedral_energy, virial);
      break;
    case math::rectangular:
      gpu::dihedral_kernel<math::rectangular><<<blocks, threads, 0, stream>>>(
          pos, dihedral_i, dihedral_j, dihedral_k, dihedral_l, dihedral_type,
          K, cospd, pd, m, atom_energy_group, num_dihedrals,
          gpu::Periodicity<math::rectangular>(box), force, dihedral_energy, virial);
      break;
    case math::triclinic:
      gpu::dihedral_kernel<math::triclinic><<<blocks, threads, 0, stream>>>(
          pos, dihedral_i, dihedral_j, dihedral_k, dihedral_l, dihedral_type,
          K, cospd, pd, m, atom_energy_group, num_dihedrals,
          gpu::Periodicity<math::triclinic>(box), force, dihedral_energy, virial);
      break;
    default:
      break;
  }
}
