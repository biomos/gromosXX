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
#include "gpu/cuda/memory/precision.h"
#include "math/gmath.h"

#include "m_shake_kernels.h"

namespace gpu {

// Fixed to match launch_m_shake_solvent's launch config below -- sized
// for the shared-memory virial reduction (9 components x THREADS
// partial sums), not a hard kernel limit.
constexpr unsigned M_SHAKE_THREADS = 128;

__global__ void m_shake_solvent_kernel(
    FPL3_TYPE* __restrict__ pos,
    const FPL3_TYPE* __restrict__ old_pos,
    const gpu::MShakeConstraint* __restrict__ constr,
    const FPL_TYPE* __restrict__ factor,
    const FPL_TYPE* __restrict__ constr_length2,
    const FPL_TYPE* __restrict__ mass_i,
    unsigned first_atom,
    unsigned num_molecules,
    FPL_TYPE tolerance,
    unsigned max_iterations,
    FPL_TYPE dt2i,
    bool do_virial,
    FPH3_TYPE* __restrict__ constraint_force,
    double* __restrict__ virial,
    int* __restrict__ error_flag) {

  const unsigned mol = blockIdx.x * blockDim.x + threadIdx.x;
  const bool active = mol < num_molecules;

  // Virial stays high precision (matches the global accumulator it
  // feeds), even though the geometry/matrix math it's built from is
  // FPL_TYPE -- same convention as shake_kernels.cu.
  double v_local[9] = {0,0,0,0,0,0,0,0,0};

  if (active) {
    const unsigned base = first_atom + mol * 3u;

    // dist_old never changes across iterations (old_pos is this step's
    // fixed reference) -- computed once, unlike dist_new below.
    FPL3_TYPE dist_old[3];
    for (unsigned k = 0; k < 3; ++k) {
      dist_old[k] = old_pos[base + constr[k].i] - old_pos[base + constr[k].j];
    }

    FPL3_TYPE cf0 = make_FPL3(FPL_TYPE(0), FPL_TYPE(0), FPL_TYPE(0));
    FPL3_TYPE cf1 = make_FPL3(FPL_TYPE(0), FPL_TYPE(0), FPL_TYPE(0));
    FPL3_TYPE cf2 = make_FPL3(FPL_TYPE(0), FPL_TYPE(0), FPL_TYPE(0));

    const FPL_TYPE tol2 = tolerance * FPL_TYPE(2);
    unsigned iterations = 0;
    bool convergence = false;

    while (!convergence) {
      convergence = true;

      FPL3_TYPE dist_new[3];
      for (unsigned k = 0; k < 3; ++k) {
        dist_new[k] = pos[base + constr[k].i] - pos[base + constr[k].j];
      }
      const FPL_TYPE dist2_0 = dot(dist_new[0], dist_new[0]);
      const FPL_TYPE dist2_1 = dot(dist_new[1], dist_new[1]);
      const FPL_TYPE dist2_2 = dot(dist_new[2], dist_new[2]);
      const FPL_TYPE diff0 = constr_length2[0] - dist2_0;
      const FPL_TYPE diff1 = constr_length2[1] - dist2_1;
      const FPL_TYPE diff2 = constr_length2[2] - dist2_2;

      if (fabs(diff0) >= constr_length2[0] * tol2 ||
          fabs(diff1) >= constr_length2[1] * tol2 ||
          fabs(diff2) >= constr_length2[2] * tol2) {

        // A(k,l) = dot(dist_old[l], dist_new[k]) * factor(k,l),
        // factor row-major: factor[3*k+l].
        FPL_TYPE A[9];
        for (unsigned k = 0; k < 3; ++k) {
          for (unsigned l = 0; l < 3; ++l) {
            A[3*k+l] = dot(dist_old[l], dist_new[k]) * factor[3*k+l];
          }
        }

        if (A[0] < constr_length2[0] * FPL_TYPE(1.0e-12) ||
            A[4] < constr_length2[1] * FPL_TYPE(1.0e-12) ||
            A[8] < constr_length2[2] * FPL_TYPE(1.0e-12)) {
          atomicExch(error_flag, 1);
          return;
        }

        // 3x3 cofactor inverse, matching math::inverse (gmath.h).
        const FPL_TYPE det = A[0]*(A[4]*A[8]-A[5]*A[7])
                          - A[1]*(A[3]*A[8]-A[5]*A[6])
                          + A[2]*(A[3]*A[7]-A[4]*A[6]);
        const FPL_TYPE idet = FPL_TYPE(1) / det;
        const FPL_TYPE Ai00 = (A[4]*A[8]-A[5]*A[7]) * idet;
        const FPL_TYPE Ai01 = (A[2]*A[7]-A[1]*A[8]) * idet;
        const FPL_TYPE Ai02 = (A[1]*A[5]-A[2]*A[4]) * idet;
        const FPL_TYPE Ai10 = (A[5]*A[6]-A[3]*A[8]) * idet;
        const FPL_TYPE Ai11 = (A[0]*A[8]-A[2]*A[6]) * idet;
        const FPL_TYPE Ai12 = (A[2]*A[3]-A[0]*A[5]) * idet;
        const FPL_TYPE Ai20 = (A[3]*A[7]-A[4]*A[6]) * idet;
        const FPL_TYPE Ai21 = (A[1]*A[6]-A[0]*A[7]) * idet;
        const FPL_TYPE Ai22 = (A[0]*A[4]-A[1]*A[3]) * idet;

        const FPL_TYPE f0 = (Ai00*diff0 + Ai01*diff1 + Ai02*diff2) * FPL_TYPE(0.5);
        const FPL_TYPE f1 = (Ai10*diff0 + Ai11*diff1 + Ai12*diff2) * FPL_TYPE(0.5);
        const FPL_TYPE f2 = (Ai20*diff0 + Ai21*diff1 + Ai22*diff2) * FPL_TYPE(0.5);

        const FPL3_TYPE f01 = f0 * dist_old[0] + f1 * dist_old[1];
        const FPL3_TYPE f02 = f2 * dist_old[2] - f0 * dist_old[0];
        const FPL3_TYPE f12 = f1 * dist_old[1] + f2 * dist_old[2];

        cf0 += f01;
        cf1 += f02;
        cf2 -= f12;

        if (do_virial) {
          const double da[3] = {static_cast<double>(dist_old[0].x),
                                 static_cast<double>(dist_old[0].y),
                                 static_cast<double>(dist_old[0].z)};
          const double db[3] = {static_cast<double>(dist_old[1].x),
                                 static_cast<double>(dist_old[1].y),
                                 static_cast<double>(dist_old[1].z)};
          const double dc[3] = {static_cast<double>(dist_old[2].x),
                                 static_cast<double>(dist_old[2].y),
                                 static_cast<double>(dist_old[2].z)};
          const double f0d = static_cast<double>(f0);
          const double f1d = static_cast<double>(f1);
          const double f2d = static_cast<double>(f2);
          const double dt2i_d = static_cast<double>(dt2i);
          for (unsigned a = 0; a < 3; ++a) {
            for (unsigned aa = 0; aa < 3; ++aa) {
              v_local[3*a+aa] -= (da[a]*da[aa]*f0d + db[a]*db[aa]*f1d + dc[a]*dc[aa]*f2d) * dt2i_d;
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

    // cf0/cf1/cf2 accumulated in FPL_TYPE (register, throughput-
    // critical); widen only at this final write into the shared FPH
    // accumulator, same convention as shake_kernels.cu.
    constraint_force[base + 0].x += static_cast<double>(cf0.x);
    constraint_force[base + 0].y += static_cast<double>(cf0.y);
    constraint_force[base + 0].z += static_cast<double>(cf0.z);
    constraint_force[base + 1].x += static_cast<double>(cf1.x);
    constraint_force[base + 1].y += static_cast<double>(cf1.y);
    constraint_force[base + 1].z += static_cast<double>(cf1.z);
    constraint_force[base + 2].x += static_cast<double>(cf2.x);
    constraint_force[base + 2].y += static_cast<double>(cf2.y);
    constraint_force[base + 2].z += static_cast<double>(cf2.z);
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
    FPL3_TYPE* pos,
    const FPL3_TYPE* old_pos,
    const gpu::MShakeConstraint* constr,
    const FPL_TYPE* factor,
    const FPL_TYPE* constr_length2,
    const FPL_TYPE* mass_i,
    unsigned first_atom,
    unsigned num_molecules,
    FPL_TYPE tolerance,
    unsigned max_iterations,
    FPL_TYPE dt2i,
    bool do_virial,
    FPH3_TYPE* constraint_force,
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
