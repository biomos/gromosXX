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
 * @file m_shake_kernels.cu
 * Per-solvent-molecule M-SHAKE kernel. See m_shake_kernels.h.
 */

#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/memory/types.h"
#include "gpu/cuda/memory/cuvector.h"
#include "math/gmath.h"

#include "m_shake_kernels.h"

namespace gpu {

// Fixed to match launch_m_shake_solvent's launch config below -- sized
// for the shared-memory virial reduction (9 components x THREADS
// partial sums), not a hard kernel limit.
constexpr unsigned M_SHAKE_THREADS = 128;

__global__ void m_shake_solvent_kernel(
    double3* __restrict__ pos,
    const double3* __restrict__ old_pos,
    const gpu::MShakeConstraint* __restrict__ constr,
    const double* __restrict__ factor,
    const double* __restrict__ constr_length2,
    const double* __restrict__ mass_i,
    unsigned first_atom,
    unsigned num_molecules,
    double tolerance,
    unsigned max_iterations,
    double dt2i,
    bool do_virial,
    double3* __restrict__ constraint_force,
    double* __restrict__ virial,
    int* __restrict__ error_flag) {

  const unsigned mol = blockIdx.x * blockDim.x + threadIdx.x;
  const bool active = mol < num_molecules;

  double v_local[9] = {0,0,0,0,0,0,0,0,0};

  if (active) {
    const unsigned base = first_atom + mol * 3u;

    // dist_old never changes across iterations (old_pos is this step's
    // fixed reference) -- computed once, unlike dist_new below.
    double3 dist_old[3];
    for (unsigned k = 0; k < 3; ++k) {
      dist_old[k] = old_pos[base + constr[k].i] - old_pos[base + constr[k].j];
    }

    double3 cf0 = make_double3(0.0, 0.0, 0.0);
    double3 cf1 = make_double3(0.0, 0.0, 0.0);
    double3 cf2 = make_double3(0.0, 0.0, 0.0);

    const double tol2 = tolerance * 2.0;
    unsigned iterations = 0;
    bool convergence = false;

    while (!convergence) {
      convergence = true;

      double3 dist_new[3];
      for (unsigned k = 0; k < 3; ++k) {
        dist_new[k] = pos[base + constr[k].i] - pos[base + constr[k].j];
      }
      const double dist2_0 = dot(dist_new[0], dist_new[0]);
      const double dist2_1 = dot(dist_new[1], dist_new[1]);
      const double dist2_2 = dot(dist_new[2], dist_new[2]);
      const double diff0 = constr_length2[0] - dist2_0;
      const double diff1 = constr_length2[1] - dist2_1;
      const double diff2 = constr_length2[2] - dist2_2;

      if (fabs(diff0) >= constr_length2[0] * tol2 ||
          fabs(diff1) >= constr_length2[1] * tol2 ||
          fabs(diff2) >= constr_length2[2] * tol2) {

        // A(k,l) = dot(dist_old[l], dist_new[k]) * factor(k,l),
        // factor row-major: factor[3*k+l].
        double A[9];
        for (unsigned k = 0; k < 3; ++k) {
          for (unsigned l = 0; l < 3; ++l) {
            A[3*k+l] = dot(dist_old[l], dist_new[k]) * factor[3*k+l];
          }
        }

        if (A[0] < constr_length2[0] * 1.0e-12 ||
            A[4] < constr_length2[1] * 1.0e-12 ||
            A[8] < constr_length2[2] * 1.0e-12) {
          atomicExch(error_flag, 1);
          return;
        }

        // 3x3 cofactor inverse, matching math::inverse (gmath.h).
        const double det = A[0]*(A[4]*A[8]-A[5]*A[7])
                          - A[1]*(A[3]*A[8]-A[5]*A[6])
                          + A[2]*(A[3]*A[7]-A[4]*A[6]);
        const double idet = 1.0 / det;
        const double Ai00 = (A[4]*A[8]-A[5]*A[7]) * idet;
        const double Ai01 = (A[2]*A[7]-A[1]*A[8]) * idet;
        const double Ai02 = (A[1]*A[5]-A[2]*A[4]) * idet;
        const double Ai10 = (A[5]*A[6]-A[3]*A[8]) * idet;
        const double Ai11 = (A[0]*A[8]-A[2]*A[6]) * idet;
        const double Ai12 = (A[2]*A[3]-A[0]*A[5]) * idet;
        const double Ai20 = (A[3]*A[7]-A[4]*A[6]) * idet;
        const double Ai21 = (A[1]*A[6]-A[0]*A[7]) * idet;
        const double Ai22 = (A[0]*A[4]-A[1]*A[3]) * idet;

        const double f0 = (Ai00*diff0 + Ai01*diff1 + Ai02*diff2) * 0.5;
        const double f1 = (Ai10*diff0 + Ai11*diff1 + Ai12*diff2) * 0.5;
        const double f2 = (Ai20*diff0 + Ai21*diff1 + Ai22*diff2) * 0.5;

        const double3 f01 = f0 * dist_old[0] + f1 * dist_old[1];
        const double3 f02 = f2 * dist_old[2] - f0 * dist_old[0];
        const double3 f12 = f1 * dist_old[1] + f2 * dist_old[2];

        cf0 += f01;
        cf1 += f02;
        cf2 -= f12;

        if (do_virial) {
          for (unsigned a = 0; a < 3; ++a) {
            const double da0 = (a==0)?dist_old[0].x:(a==1)?dist_old[0].y:dist_old[0].z;
            const double da1 = (a==0)?dist_old[1].x:(a==1)?dist_old[1].y:dist_old[1].z;
            const double da2 = (a==0)?dist_old[2].x:(a==1)?dist_old[2].y:dist_old[2].z;
            for (unsigned aa = 0; aa < 3; ++aa) {
              const double db0 = (aa==0)?dist_old[0].x:(aa==1)?dist_old[0].y:dist_old[0].z;
              const double db1 = (aa==0)?dist_old[1].x:(aa==1)?dist_old[1].y:dist_old[1].z;
              const double db2 = (aa==0)?dist_old[2].x:(aa==1)?dist_old[2].y:dist_old[2].z;
              v_local[3*a+aa] -= (da0*db0*f0 + da1*db1*f1 + da2*db2*f2) * dt2i;
            }
          }
        }

        pos[base + 0] += f01 * mass_i[0];
        pos[base + 1] += f02 * mass_i[1];
        pos[base + 2] -= f12 * mass_i[2];

        convergence = false;
      }

      ++iterations;
      if (iterations > max_iterations) {
        atomicExch(error_flag, 2);
        return;
      }
    }

    constraint_force[base + 0] += cf0;
    constraint_force[base + 1] += cf1;
    constraint_force[base + 2] += cf2;
  }

  if (do_virial) {
    __shared__ double s_v[9][M_SHAKE_THREADS];
    for (unsigned k = 0; k < 9; ++k) s_v[k][threadIdx.x] = v_local[k];
    __syncthreads();
    if (threadIdx.x == 0) {
      for (unsigned k = 0; k < 9; ++k) {
        double sum = 0.0;
        for (unsigned t = 0; t < blockDim.x; ++t) sum += s_v[k][t];
        if (sum != 0.0) atomicAdd(&virial[k], sum);
      }
    }
  }
}

} // namespace gpu

void gpu::launch_m_shake_solvent(
    double3* pos,
    const double3* old_pos,
    const gpu::MShakeConstraint* constr,
    const double* factor,
    const double* constr_length2,
    const double* mass_i,
    unsigned first_atom,
    unsigned num_molecules,
    double tolerance,
    unsigned max_iterations,
    double dt2i,
    bool do_virial,
    double3* constraint_force,
    double* virial,
    int* error_flag,
    cudaStream_t stream) {

  if (num_molecules == 0) return;

  const unsigned threads = gpu::M_SHAKE_THREADS;
  const unsigned blocks = (num_molecules + threads - 1) / threads;

  gpu::m_shake_solvent_kernel<<<blocks, threads, 0, stream>>>(
      pos, old_pos, constr, factor, constr_length2, mass_i,
      first_atom, num_molecules, tolerance, max_iterations, dt2i, do_virial,
      constraint_force, virial, error_flag);
}
