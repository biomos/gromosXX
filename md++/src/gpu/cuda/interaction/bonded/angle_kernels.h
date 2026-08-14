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
 * @file angle_kernels.h
 * Host-callable entry point for the harmonic (cosine) bond-angle
 * force/energy kernel (PLAN.md §10 step 14, bonded forces). One thread
 * per angle term, direct global atomicAdd -- same shape and rationale
 * as quartic_bond_kernels.h (small, static term list; no pairlist).
 *
 * Exact CPU formula (angle_interaction.cc's
 * _calculate_angle_interactions): for angle term (i, j, k, type),
 * rij = nearest_image(pos(i), pos(j)), rkj = nearest_image(pos(k), pos(j)),
 * dij = |rij|, dkj = |rkj|, cost = dot(rij,rkj)/(dij*dkj),
 * df = -K[type] * (cost - cos0[type]),
 * fi = df/dij * (rkj/dkj - rij/dij*cost), fk = df/dkj * (rij/dij - rkj/dkj*cost),
 * fj = -fi - fk, force(i)+=fi, force(j)+=fj, force(k)+=fk,
 * e = 0.5*K[type]*(cost-cos0[type])^2, accumulated into
 * angle_energy[atom_energy_group[i]]. Virial is unconditional (the CPU
 * code's `if (V == math::atomic_virial)` gate is dead/commented-out):
 * virial(a,b) += rij(a)*fi(b) + rkj(a)*fk(b).
 */

#pragma once

#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/cuvector.h"
#include "math/gmath.h"

namespace gpu {

  void launch_angle(
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
      cudaStream_t stream = 0);

} // namespace gpu
