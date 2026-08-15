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
 * @file settle_kernels.cu
 * SETTLE kernel. See settle_kernels.h. Translates algorithm::Settle::
 * solvent() (settle.cc) line-for-line, one thread per molecule.
 */

#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/memory/types.h"
#include "gpu/cuda/memory/cuvector.h"
#include "math/gmath.h"

#include "settle_kernels.h"

namespace gpu {

__global__ void settle_kernel(
    double3* __restrict__ pos,
    const double3* __restrict__ old_pos,
    double3* __restrict__ vel,
    unsigned first_atom,
    unsigned num_molecules,
    double ra,
    double rb,
    double rc,
    double mass_O,
    double mass_H,
    double dt_i,
    double dt2_i,
    bool do_velocity,
    double3* __restrict__ constraint_force,
    double* __restrict__ virial,
    int* __restrict__ error_flag) {

  const unsigned mol = blockIdx.x * blockDim.x + threadIdx.x;
  if (mol >= num_molecules) return;

  const unsigned base = first_atom + mol * 3;

  const double3 pos_new0 = pos[base];
  const double3 pos_new1 = pos[base + 1];
  const double3 pos_new2 = pos[base + 2];
  const double3 pos_old0 = old_pos[base];
  const double3 pos_old1 = old_pos[base + 1];
  const double3 pos_old2 = old_pos[base + 2];

  // vectors in the plane of the old positions
  double3 b0 = pos_old1 - pos_old0;
  double3 c0 = pos_old2 - pos_old0;

  // centre of mass of new positions
  const double3 d0 = (pos_new0 * mass_O + pos_new1 * mass_H + pos_new2 * mass_H) /
                      (mass_O + mass_H + mass_H);

  // move the origin to the centre of mass
  double3 a1 = pos_new0 - d0;
  double3 b1 = pos_new1 - d0;
  double3 c1 = pos_new2 - d0;

  // vectors describing transformation from original coordinate system to
  // the centre of mass originated coordinate system
  double3 n0 = cross(b0, c0);
  double3 n1 = cross(a1, n0);
  double3 n2 = cross(n0, n1);
  n0 = n0 / sqrt(dot(n0, n0)); // this can give a NaN but it is very unlikely.
  n1 = n1 / sqrt(dot(n1, n1));
  n2 = n2 / sqrt(dot(n2, n2));

  // generate new normal vectors from the transpose in order to undo
  // the transformation
  const double3 m1{n1.x, n2.x, n0.x};
  const double3 m2{n1.y, n2.y, n0.y};
  const double3 m0{n1.z, n2.z, n0.z};

  // do the transformation to the centre of mass originated coordinate system
  // of the old positions
  b0 = double3{dot(n1, b0), dot(n2, b0), dot(n0, b0)};
  c0 = double3{dot(n1, c0), dot(n2, c0), dot(n0, c0)};

  // and of the new positions
  a1 = double3{dot(n1, a1), dot(n2, a1), dot(n0, a1)};
  b1 = double3{dot(n1, b1), dot(n2, b1), dot(n0, b1)};
  c1 = double3{dot(n1, c1), dot(n2, c1), dot(n0, c1)};

  // now we can compute positions of canonical water
  const double sinphi = a1.z / ra; // (A8)
  const double one_minus_sinphi2 = 1.0 - sinphi * sinphi;
  if (one_minus_sinphi2 < 0.0) {
    atomicExch(error_flag, 1);
    return;
  }
  const double cosphi = sqrt(one_minus_sinphi2);

  const double sinpsi = (b1.z - c1.z) / (2.0 * rc * cosphi); // (A9)
  const double one_minus_sinpsi2 = 1.0 - sinpsi * sinpsi;
  if (one_minus_sinpsi2 < 0.0) {
    atomicExch(error_flag, 1);
    return;
  }
  const double cospsi = sqrt(one_minus_sinpsi2);

  const double minus_rb_cosphi = -rb * cosphi;
  const double rc_cospsi = rc * cospsi;
  const double rc_sinpsi_sinphi = rc * sinpsi * sinphi;
  const double rc_sinpsi_cosphi = rc * sinpsi * cosphi;

  const double x_a2 = 0.0;
  const double x_b2 = -rc_cospsi;
  const double x_c2 = rc_cospsi;

  const double y_a2 = ra * cosphi;
  const double y_b2 = minus_rb_cosphi - rc_sinpsi_sinphi;
  const double y_c2 = minus_rb_cosphi + rc_sinpsi_sinphi;

  const double z_a2 = ra * sinphi; // (A5)
  const double z_b2 = -rb * sinphi + rc_sinpsi_cosphi; // (A6)
  const double z_c2 = -rb * sinphi - rc_sinpsi_cosphi; // (A7)

  const double3 a2{x_a2, y_a2, z_a2};
  const double3 b2{x_b2, y_b2, z_b2};
  const double3 c2{x_c2, y_c2, z_c2};

  const double alpha = b2.x * (b0.x - c0.x) + b0.y * b2.y + c0.y * c2.y;
  const double beta = b2.x * (c0.y - b0.y) + b0.x * b2.y + c0.x * c2.y;
  const double gamma = b0.x * b1.y - b1.x * b0.y + c0.x * c1.y - c1.x * c0.y;

  const double alpha2_beta2 = alpha * alpha + beta * beta;
  const double under_sqrt = alpha2_beta2 - gamma * gamma;
  if (under_sqrt < 0.0) {
    atomicExch(error_flag, 1);
    return;
  }
  const double sintheta = (alpha * gamma - beta * sqrt(under_sqrt)) / alpha2_beta2; // (A17)
  const double one_minus_sintheta2 = 1.0 - sintheta * sintheta;
  if (one_minus_sintheta2 < 0.0) {
    atomicExch(error_flag, 1);
    return;
  }
  const double costheta = sqrt(one_minus_sintheta2);

  const double3 a3{-a2.y * sintheta, a2.y * costheta, a1.z};
  const double3 b3{b2.x * costheta - b2.y * sintheta,
                    b2.x * sintheta + b2.y * costheta, b1.z};
  const double3 c3{-b2.x * costheta - c2.y * sintheta,
                    -b2.x * sintheta + c2.y * costheta, c1.z};

  const double3 pos_a = double3{dot(a3, m1), dot(a3, m2), dot(a3, m0)} + d0;
  const double3 pos_b = double3{dot(b3, m1), dot(b3, m2), dot(b3, m0)} + d0;
  const double3 pos_c = double3{dot(c3, m1), dot(c3, m2), dot(c3, m0)} + d0;

  const double3 d_a = pos_a - pos_new0;
  const double3 d_b = pos_b - pos_new1;
  const double3 d_c = pos_c - pos_new2;

  pos[base] = pos_a;
  pos[base + 1] = pos_b;
  pos[base + 2] = pos_c;

  const double3 cf0 = d_a * mass_O * dt2_i;
  const double3 cf1 = d_b * mass_H * dt2_i;
  const double3 cf2 = d_c * mass_H * dt2_i;
  constraint_force[base] = cf0;
  constraint_force[base + 1] = cf1;
  constraint_force[base + 2] = cf2;

  if (do_velocity) {
    vel[base] += d_a * dt_i;
    vel[base + 1] += d_b * dt_i;
    vel[base + 2] += d_c * dt_i;
  }

  // Unconditional atomic virial -- whether it's used is the host's
  // decision (matches the CPU's live pcouple.virial gate, see this
  // file's header doc comment).
  atomicAdd(&virial[0], pos_old0.x*cf0.x + pos_old1.x*cf1.x + pos_old2.x*cf2.x);
  atomicAdd(&virial[1], pos_old0.x*cf0.y + pos_old1.x*cf1.y + pos_old2.x*cf2.y);
  atomicAdd(&virial[2], pos_old0.x*cf0.z + pos_old1.x*cf1.z + pos_old2.x*cf2.z);
  atomicAdd(&virial[3], pos_old0.y*cf0.x + pos_old1.y*cf1.x + pos_old2.y*cf2.x);
  atomicAdd(&virial[4], pos_old0.y*cf0.y + pos_old1.y*cf1.y + pos_old2.y*cf2.y);
  atomicAdd(&virial[5], pos_old0.y*cf0.z + pos_old1.y*cf1.z + pos_old2.y*cf2.z);
  atomicAdd(&virial[6], pos_old0.z*cf0.x + pos_old1.z*cf1.x + pos_old2.z*cf2.x);
  atomicAdd(&virial[7], pos_old0.z*cf0.y + pos_old1.z*cf1.y + pos_old2.z*cf2.y);
  atomicAdd(&virial[8], pos_old0.z*cf0.z + pos_old1.z*cf1.z + pos_old2.z*cf2.z);
}

} // namespace gpu

void gpu::launch_settle(
    double3* pos,
    const double3* old_pos,
    double3* vel,
    unsigned first_atom,
    unsigned num_molecules,
    double mass_O,
    double mass_H,
    double dist_OH,
    double dist_HH,
    double dt_i,
    bool do_velocity,
    double3* constraint_force,
    double* virial,
    int* error_flag,
    cudaStream_t stream) {

  if (num_molecules == 0) return;

  const double half_mO_div_mH = 0.5 * mass_O / mass_H;
  const double rc = 0.5 * dist_HH;
  const double ra = sqrt(dist_OH * dist_OH - rc * rc) / (1.0 + half_mO_div_mH);
  const double rb = half_mO_div_mH * ra;
  const double dt2_i = dt_i * dt_i;

  const unsigned threads = 128;
  const unsigned blocks = (num_molecules + threads - 1) / threads;

  gpu::settle_kernel<<<blocks, threads, 0, stream>>>(
      pos, old_pos, vel, first_atom, num_molecules, ra, rb, rc,
      mass_O, mass_H, dt_i, dt2_i, do_velocity, constraint_force, virial, error_flag);
}
