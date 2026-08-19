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
 * @file lincs_kernels.cu
 * LINCS kernels. See lincs_kernels.h.
 */

#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/memory/types.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/math/periodicity.h"
#include "math/gmath.h"

#include "lincs_kernels.h"

namespace gpu {

template <math::boundary_enum BOUNDARY>
__global__ void lincs_compute_b_kernel(
    const FPL3_TYPE* __restrict__ old_pos,
    const gpu::LincsConstraint* __restrict__ constraints,
    unsigned num_constr_per_instance,
    unsigned num_instances,
    unsigned first_atom,
    unsigned atom_stride_per_instance,
    gpu::Periodicity<BOUNDARY> periodicity,
    FPL3_TYPE* __restrict__ B) {

  const unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
  const unsigned total = num_instances * num_constr_per_instance;
  if (idx >= total) return;

  const unsigned instance = idx / num_constr_per_instance;
  const unsigned local_c = idx % num_constr_per_instance;
  const unsigned base = first_atom + instance * atom_stride_per_instance;
  const unsigned ai = base + constraints[local_c].i;
  const unsigned aj = base + constraints[local_c].j;

  const FPL3_TYPE ref_r = periodicity.nearest_image(old_pos[ai], old_pos[aj]);
  const FPL_TYPE norm = sqrt(dot(ref_r, ref_r));
  B[idx] = ref_r / norm;
}

template <math::boundary_enum BOUNDARY>
__global__ void lincs_init_rhs_kernel(
    const FPL3_TYPE* __restrict__ pos,
    const gpu::LincsConstraint* __restrict__ constraints,
    unsigned num_constr_per_instance,
    unsigned num_instances,
    unsigned first_atom,
    unsigned atom_stride_per_instance,
    gpu::Periodicity<BOUNDARY> periodicity,
    const FPL3_TYPE* __restrict__ B,
    FPL_TYPE* __restrict__ rhs,
    FPL_TYPE* __restrict__ sol) {

  const unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
  const unsigned total = num_instances * num_constr_per_instance;
  if (idx >= total) return;

  const unsigned instance = idx / num_constr_per_instance;
  const unsigned local_c = idx % num_constr_per_instance;
  const unsigned base = first_atom + instance * atom_stride_per_instance;
  const unsigned ai = base + constraints[local_c].i;
  const unsigned aj = base + constraints[local_c].j;

  const FPL3_TYPE r = periodicity.nearest_image(pos[ai], pos[aj]);
  const FPL_TYPE val = constraints[local_c].sdiag * (dot(B[idx], r) - constraints[local_c].r0);
  rhs[idx] = val;
  sol[idx] = val;
}

template <math::boundary_enum BOUNDARY>
__global__ void lincs_rotation_rhs_kernel(
    const FPL3_TYPE* __restrict__ pos,
    const gpu::LincsConstraint* __restrict__ constraints,
    unsigned num_constr_per_instance,
    unsigned num_instances,
    unsigned first_atom,
    unsigned atom_stride_per_instance,
    gpu::Periodicity<BOUNDARY> periodicity,
    FPL_TYPE* __restrict__ rhs,
    FPL_TYPE* __restrict__ sol,
    int* __restrict__ rotation_count) {

  const unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
  const unsigned total = num_instances * num_constr_per_instance;
  if (idx >= total) return;

  const unsigned instance = idx / num_constr_per_instance;
  const unsigned local_c = idx % num_constr_per_instance;
  const unsigned base = first_atom + instance * atom_stride_per_instance;
  const unsigned ai = base + constraints[local_c].i;
  const unsigned aj = base + constraints[local_c].j;

  const FPL3_TYPE r = periodicity.nearest_image(pos[ai], pos[aj]);
  const FPL_TYPE r0 = constraints[local_c].r0;
  const FPL_TYPE diff = FPL_TYPE(2) * r0 * r0 - dot(r, r);
  FPL_TYPE p = FPL_TYPE(0);
  if (diff > FPL_TYPE(0)) {
    p = sqrt(diff);
  } else {
    atomicAdd(rotation_count, 1);
  }
  const FPL_TYPE val = constraints[local_c].sdiag * (r0 - p);
  rhs[idx] = val;
  sol[idx] = val;
}

__global__ void lincs_round_kernel(
    const FPL3_TYPE* __restrict__ B,
    const unsigned* __restrict__ coupled_offset,
    const unsigned* __restrict__ coupled_index,
    const FPL_TYPE* __restrict__ coupled_coef,
    unsigned num_constr_per_instance,
    unsigned num_instances,
    const FPL_TYPE* __restrict__ rhs_in,
    FPL_TYPE* __restrict__ rhs_out,
    FPL_TYPE* __restrict__ sol) {

  const unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
  const unsigned total = num_instances * num_constr_per_instance;
  if (idx >= total) return;

  const unsigned instance = idx / num_constr_per_instance;
  const unsigned local_c = idx % num_constr_per_instance;
  const FPL3_TYPE Bi = B[idx];

  FPL_TYPE acc = FPL_TYPE(0);
  const unsigned start = coupled_offset[local_c];
  const unsigned end = coupled_offset[local_c + 1];
  for (unsigned k = start; k < end; ++k) {
    const unsigned global_coupled = instance * num_constr_per_instance + coupled_index[k];
    acc += coupled_coef[k] * dot(Bi, B[global_coupled]) * rhs_in[global_coupled];
  }
  rhs_out[idx] = acc;
  sol[idx] += acc;
}

__global__ void lincs_apply_kernel(
    FPL3_TYPE* __restrict__ pos,
    const gpu::LincsConstraint* __restrict__ constraints,
    const FPL3_TYPE* __restrict__ B,
    const FPL_TYPE* __restrict__ sol,
    unsigned num_constr_per_instance,
    unsigned num_instances,
    unsigned first_atom,
    unsigned atom_stride_per_instance) {

  const unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
  const unsigned total = num_instances * num_constr_per_instance;
  if (idx >= total) return;

  const unsigned instance = idx / num_constr_per_instance;
  const unsigned local_c = idx % num_constr_per_instance;
  const unsigned base = first_atom + instance * atom_stride_per_instance;
  const unsigned ai = base + constraints[local_c].i;
  const unsigned aj = base + constraints[local_c].j;

  const FPL3_TYPE Bi = B[idx];
  const FPL_TYPE coeff = constraints[local_c].sdiag * sol[idx];
  const FPL3_TYPE corr_i = Bi / constraints[local_c].mass_i * coeff;
  const FPL3_TYPE corr_j = Bi / constraints[local_c].mass_j * coeff;

  atomicAdd(&pos[ai].x, -corr_i.x);
  atomicAdd(&pos[ai].y, -corr_i.y);
  atomicAdd(&pos[ai].z, -corr_i.z);
  atomicAdd(&pos[aj].x, corr_j.x);
  atomicAdd(&pos[aj].y, corr_j.y);
  atomicAdd(&pos[aj].z, corr_j.z);
}

} // namespace gpu

namespace {
  inline unsigned num_blocks(unsigned total, unsigned threads) {
    return (total + threads - 1) / threads;
  }
  constexpr unsigned kThreads = 128;
}

void gpu::launch_lincs_compute_b(
    const FPL3_TYPE* old_pos,
    const gpu::LincsConstraint* constraints,
    unsigned num_constr_per_instance,
    unsigned num_instances,
    unsigned first_atom,
    unsigned atom_stride_per_instance,
    math::boundary_enum boundary,
    math::Box box,
    FPL3_TYPE* B,
    cudaStream_t stream) {

  const unsigned total = num_instances * num_constr_per_instance;
  if (total == 0) return;
  const unsigned blocks = num_blocks(total, kThreads);

  switch (boundary) {
    case math::vacuum:
      gpu::lincs_compute_b_kernel<math::vacuum><<<blocks, kThreads, 0, stream>>>(
          old_pos, constraints, num_constr_per_instance, num_instances, first_atom,
          atom_stride_per_instance, gpu::Periodicity<math::vacuum>(box), B);
      break;
    case math::rectangular:
      gpu::lincs_compute_b_kernel<math::rectangular><<<blocks, kThreads, 0, stream>>>(
          old_pos, constraints, num_constr_per_instance, num_instances, first_atom,
          atom_stride_per_instance, gpu::Periodicity<math::rectangular>(box), B);
      break;
    case math::triclinic:
      gpu::lincs_compute_b_kernel<math::triclinic><<<blocks, kThreads, 0, stream>>>(
          old_pos, constraints, num_constr_per_instance, num_instances, first_atom,
          atom_stride_per_instance, gpu::Periodicity<math::triclinic>(box), B);
      break;
    default:
      break;
  }
}

void gpu::launch_lincs_init_rhs(
    const FPL3_TYPE* pos,
    const gpu::LincsConstraint* constraints,
    unsigned num_constr_per_instance,
    unsigned num_instances,
    unsigned first_atom,
    unsigned atom_stride_per_instance,
    math::boundary_enum boundary,
    math::Box box,
    const FPL3_TYPE* B,
    FPL_TYPE* rhs,
    FPL_TYPE* sol,
    cudaStream_t stream) {

  const unsigned total = num_instances * num_constr_per_instance;
  if (total == 0) return;
  const unsigned blocks = num_blocks(total, kThreads);

  switch (boundary) {
    case math::vacuum:
      gpu::lincs_init_rhs_kernel<math::vacuum><<<blocks, kThreads, 0, stream>>>(
          pos, constraints, num_constr_per_instance, num_instances, first_atom,
          atom_stride_per_instance, gpu::Periodicity<math::vacuum>(box), B, rhs, sol);
      break;
    case math::rectangular:
      gpu::lincs_init_rhs_kernel<math::rectangular><<<blocks, kThreads, 0, stream>>>(
          pos, constraints, num_constr_per_instance, num_instances, first_atom,
          atom_stride_per_instance, gpu::Periodicity<math::rectangular>(box), B, rhs, sol);
      break;
    case math::triclinic:
      gpu::lincs_init_rhs_kernel<math::triclinic><<<blocks, kThreads, 0, stream>>>(
          pos, constraints, num_constr_per_instance, num_instances, first_atom,
          atom_stride_per_instance, gpu::Periodicity<math::triclinic>(box), B, rhs, sol);
      break;
    default:
      break;
  }
}

void gpu::launch_lincs_rotation_rhs(
    const FPL3_TYPE* pos,
    const gpu::LincsConstraint* constraints,
    unsigned num_constr_per_instance,
    unsigned num_instances,
    unsigned first_atom,
    unsigned atom_stride_per_instance,
    math::boundary_enum boundary,
    math::Box box,
    FPL_TYPE* rhs,
    FPL_TYPE* sol,
    int* rotation_count,
    cudaStream_t stream) {

  const unsigned total = num_instances * num_constr_per_instance;
  if (total == 0) return;
  const unsigned blocks = num_blocks(total, kThreads);

  switch (boundary) {
    case math::vacuum:
      gpu::lincs_rotation_rhs_kernel<math::vacuum><<<blocks, kThreads, 0, stream>>>(
          pos, constraints, num_constr_per_instance, num_instances, first_atom,
          atom_stride_per_instance, gpu::Periodicity<math::vacuum>(box), rhs, sol, rotation_count);
      break;
    case math::rectangular:
      gpu::lincs_rotation_rhs_kernel<math::rectangular><<<blocks, kThreads, 0, stream>>>(
          pos, constraints, num_constr_per_instance, num_instances, first_atom,
          atom_stride_per_instance, gpu::Periodicity<math::rectangular>(box), rhs, sol, rotation_count);
      break;
    case math::triclinic:
      gpu::lincs_rotation_rhs_kernel<math::triclinic><<<blocks, kThreads, 0, stream>>>(
          pos, constraints, num_constr_per_instance, num_instances, first_atom,
          atom_stride_per_instance, gpu::Periodicity<math::triclinic>(box), rhs, sol, rotation_count);
      break;
    default:
      break;
  }
}

void gpu::launch_lincs_round(
    const FPL3_TYPE* B,
    const unsigned* coupled_offset,
    const unsigned* coupled_index,
    const FPL_TYPE* coupled_coef,
    unsigned num_constr_per_instance,
    unsigned num_instances,
    const FPL_TYPE* rhs_in,
    FPL_TYPE* rhs_out,
    FPL_TYPE* sol,
    cudaStream_t stream) {

  const unsigned total = num_instances * num_constr_per_instance;
  if (total == 0) return;
  const unsigned blocks = num_blocks(total, kThreads);
  gpu::lincs_round_kernel<<<blocks, kThreads, 0, stream>>>(
      B, coupled_offset, coupled_index, coupled_coef, num_constr_per_instance,
      num_instances, rhs_in, rhs_out, sol);
}

void gpu::launch_lincs_apply(
    FPL3_TYPE* pos,
    const gpu::LincsConstraint* constraints,
    const FPL3_TYPE* B,
    const FPL_TYPE* sol,
    unsigned num_constr_per_instance,
    unsigned num_instances,
    unsigned first_atom,
    unsigned atom_stride_per_instance,
    cudaStream_t stream) {

  const unsigned total = num_instances * num_constr_per_instance;
  if (total == 0) return;
  const unsigned blocks = num_blocks(total, kThreads);
  gpu::lincs_apply_kernel<<<blocks, kThreads, 0, stream>>>(
      pos, constraints, B, sol, num_constr_per_instance, num_instances,
      first_atom, atom_stride_per_instance);
}
