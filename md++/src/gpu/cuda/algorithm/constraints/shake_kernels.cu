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
 * @file shake_kernels.cu
 * Per-solvent-molecule SHAKE kernel. See shake_kernels.h.
 */

#include "gpu/cuda/cuheader.h"

#include "gpu/cuda/memory/types.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/math/periodicity.h"
#include "math/gmath.h"

#include "shake_kernels.h"

namespace gpu {

template <math::boundary_enum BOUNDARY>
__global__ void shake_solvent_kernel(
    double3* __restrict__ pos,
    const double3* __restrict__ old_pos,
    const gpu::ShakeConstraint* __restrict__ constraints,
    unsigned num_constraints,
    const double* __restrict__ inv_mass_local,
    unsigned num_atoms_per_molecule,
    unsigned first_atom,
    unsigned num_molecules,
    double tolerance,
    unsigned max_iterations,
    gpu::Periodicity<BOUNDARY> periodicity,
    double dt2,
    double3* __restrict__ constraint_force,
    double* __restrict__ virial,
    int* __restrict__ error_flag) {

  const unsigned mol = blockIdx.x * blockDim.x + threadIdx.x;
  if (mol >= num_molecules) return;

  const unsigned base = first_atom + mol * num_atoms_per_molecule;

  // Cache this molecule's atoms in registers/local memory for the whole
  // Jacobi loop instead of re-reading pos[]/old_pos[] from global memory
  // on every constraint check of every iteration (measured ~12 iterations
  // per molecule -- each one previously re-fetched every atom involved).
  // old_pos never changes within apply(), so it's loaded once and never
  // written back; pos and the constraint-force accumulator are written
  // back exactly once, after convergence, instead of per-update.
  double3 local_pos[gpu::MAX_SHAKE_ATOMS_PER_MOLECULE];
  double3 local_old_pos[gpu::MAX_SHAKE_ATOMS_PER_MOLECULE];
  double3 local_cf[gpu::MAX_SHAKE_ATOMS_PER_MOLECULE];
  for (unsigned a = 0; a < num_atoms_per_molecule; ++a) {
    local_pos[a] = pos[base + a];
    local_old_pos[a] = old_pos[base + a];
    local_cf[a] = make_double3(0.0, 0.0, 0.0);
  }

  bool skip_now[gpu::MAX_SHAKE_ATOMS_PER_MOLECULE];
  bool skip_next[gpu::MAX_SHAKE_ATOMS_PER_MOLECULE];
  for (unsigned a = 0; a < num_atoms_per_molecule; ++a) {
    skip_now[a] = false;
    skip_next[a] = true;
  }

  unsigned iterations = 0;
  bool convergence = false;
  double v_local[9] = {0,0,0,0,0,0,0,0,0};

  while (!convergence) {
    convergence = true;

    for (unsigned c = 0; c < num_constraints; ++c) {
      const unsigned li = constraints[c].i;
      const unsigned lj = constraints[c].j;

      if (skip_now[li] && skip_now[lj]) continue;
      if (inv_mass_local[li] == 0.0 && inv_mass_local[lj] == 0.0) continue;

      const double3 r = periodicity.nearest_image(local_pos[li], local_pos[lj]);
      const double dist2 = dot(r, r);
      const double r0sq = constraints[c].r0sq;
      const double diff = r0sq - dist2;

      if (fabs(diff) >= r0sq * tolerance * 2.0) {
        double3 ref_r = periodicity.nearest_image(local_old_pos[li], local_old_pos[lj]);
        const double sp = dot(ref_r, r);

        // 1e-12, matching math::epsilon (math/gmath.h) -- that constant
        // isn't device-accessible (a plain non-constexpr host global),
        // so the literal is inlined here instead.
        if (sp < r0sq * 1.0e-12) {
          atomicExch(error_flag, 1);
          return;
        }

        const double lambda = diff / (sp * 2.0 * (inv_mass_local[li] + inv_mass_local[lj]));

        const double3 cons_force = lambda * ref_r;
        local_cf[li] += cons_force;
        local_cf[lj] -= cons_force;

        const double inv_dt2 = lambda / dt2;
        v_local[0] += ref_r.x * ref_r.x * inv_dt2; // (0,0)
        v_local[1] += ref_r.x * ref_r.y * inv_dt2; // (0,1)
        v_local[2] += ref_r.x * ref_r.z * inv_dt2; // (0,2)
        v_local[3] += ref_r.y * ref_r.x * inv_dt2; // (1,0)
        v_local[4] += ref_r.y * ref_r.y * inv_dt2; // (1,1)
        v_local[5] += ref_r.y * ref_r.z * inv_dt2; // (1,2)
        v_local[6] += ref_r.z * ref_r.x * inv_dt2; // (2,0)
        v_local[7] += ref_r.z * ref_r.y * inv_dt2; // (2,1)
        v_local[8] += ref_r.z * ref_r.z * inv_dt2; // (2,2)

        ref_r *= lambda;
        local_pos[li] += ref_r * inv_mass_local[li];
        local_pos[lj] -= ref_r * inv_mass_local[lj];

        convergence = false;
        skip_next[li] = false;
        skip_next[lj] = false;
      }
    }

    ++iterations;
    if (iterations > max_iterations) {
      atomicExch(error_flag, 2);
      return;
    }

    for (unsigned a = 0; a < num_atoms_per_molecule; ++a) {
      skip_now[a] = skip_next[a];
      skip_next[a] = true;
    }
  }

  for (unsigned a = 0; a < num_atoms_per_molecule; ++a) {
    pos[base + a] = local_pos[a];
    constraint_force[base + a] += local_cf[a];
  }

  // Accumulated locally across every constraint/iteration this thread
  // touched, then flushed with one atomicAdd per component instead of
  // one per constraint update -- thousands of molecule-threads otherwise
  // hammer the same 9 global addresses every iteration, serializing the
  // whole kernel (measured: ~150ms/call before this change).
  for (unsigned k = 0; k < 9; ++k) {
    if (v_local[k] != 0.0) atomicAdd(&virial[k], v_local[k]);
  }
}

} // namespace gpu

void gpu::launch_shake_solvent(
    double3* pos,
    const double3* old_pos,
    const gpu::ShakeConstraint* constraints,
    unsigned num_constraints,
    const double* inv_mass_local,
    unsigned num_atoms_per_molecule,
    unsigned first_atom,
    unsigned num_molecules,
    double tolerance,
    unsigned max_iterations,
    math::boundary_enum boundary,
    math::Box box,
    double dt2,
    double3* constraint_force,
    double* virial,
    int* error_flag,
    cudaStream_t stream) {

  if (num_molecules == 0) return;

  const unsigned threads = 128;
  const unsigned blocks = (num_molecules + threads - 1) / threads;

  switch (boundary) {
    case math::vacuum:
      gpu::shake_solvent_kernel<math::vacuum><<<blocks, threads, 0, stream>>>(
          pos, old_pos, constraints, num_constraints, inv_mass_local,
          num_atoms_per_molecule, first_atom, num_molecules, tolerance,
          max_iterations, gpu::Periodicity<math::vacuum>(box), dt2,
          constraint_force, virial, error_flag);
      break;
    case math::rectangular:
      gpu::shake_solvent_kernel<math::rectangular><<<blocks, threads, 0, stream>>>(
          pos, old_pos, constraints, num_constraints, inv_mass_local,
          num_atoms_per_molecule, first_atom, num_molecules, tolerance,
          max_iterations, gpu::Periodicity<math::rectangular>(box), dt2,
          constraint_force, virial, error_flag);
      break;
    case math::triclinic:
      gpu::shake_solvent_kernel<math::triclinic><<<blocks, threads, 0, stream>>>(
          pos, old_pos, constraints, num_constraints, inv_mass_local,
          num_atoms_per_molecule, first_atom, num_molecules, tolerance,
          max_iterations, gpu::Periodicity<math::triclinic>(box), dt2,
          constraint_force, virial, error_flag);
      break;
    default:
      break;
  }
}

namespace gpu {

template <math::boundary_enum BOUNDARY>
__global__ void shake_solute_round_kernel(
    const double3* __restrict__ pos,
    const double3* __restrict__ old_pos,
    const gpu::ShakeConstraint* __restrict__ constraints,
    unsigned num_constraints,
    const double* __restrict__ inv_mass,
    double tolerance,
    gpu::Periodicity<BOUNDARY> periodicity,
    double dt2,
    double3* __restrict__ delta,
    double3* __restrict__ constraint_force,
    double* __restrict__ virial,
    int* __restrict__ changed_flag,
    int* __restrict__ error_flag) {

  const unsigned c = blockIdx.x * blockDim.x + threadIdx.x;
  if (c >= num_constraints) return;

  const unsigned i = constraints[c].i;
  const unsigned j = constraints[c].j;

  const double3 r = periodicity.nearest_image(pos[i], pos[j]);
  const double dist2 = dot(r, r);
  const double r0sq = constraints[c].r0sq;
  const double diff = r0sq - dist2;

  if (fabs(diff) < r0sq * tolerance * 2.0) return;

  double3 ref_r = periodicity.nearest_image(old_pos[i], old_pos[j]);
  const double sp = dot(ref_r, r);

  // 1e-12, matching math::epsilon -- see shake_solvent_kernel's comment.
  if (sp < r0sq * 1.0e-12) {
    atomicExch(error_flag, 1);
    return;
  }

  const double lambda = diff / (sp * 2.0 * (inv_mass[i] + inv_mass[j]));

  const double3 cons_force = lambda * ref_r;
  atomicAdd(&constraint_force[i].x, cons_force.x);
  atomicAdd(&constraint_force[i].y, cons_force.y);
  atomicAdd(&constraint_force[i].z, cons_force.z);
  atomicAdd(&constraint_force[j].x, -cons_force.x);
  atomicAdd(&constraint_force[j].y, -cons_force.y);
  atomicAdd(&constraint_force[j].z, -cons_force.z);

  atomicAdd(&virial[0], ref_r.x * ref_r.x * lambda / dt2);
  atomicAdd(&virial[1], ref_r.x * ref_r.y * lambda / dt2);
  atomicAdd(&virial[2], ref_r.x * ref_r.z * lambda / dt2);
  atomicAdd(&virial[3], ref_r.y * ref_r.x * lambda / dt2);
  atomicAdd(&virial[4], ref_r.y * ref_r.y * lambda / dt2);
  atomicAdd(&virial[5], ref_r.y * ref_r.z * lambda / dt2);
  atomicAdd(&virial[6], ref_r.z * ref_r.x * lambda / dt2);
  atomicAdd(&virial[7], ref_r.z * ref_r.y * lambda / dt2);
  atomicAdd(&virial[8], ref_r.z * ref_r.z * lambda / dt2);

  ref_r *= lambda;
  const double3 di = ref_r * inv_mass[i];
  const double3 dj = ref_r * inv_mass[j];
  atomicAdd(&delta[i].x, di.x);
  atomicAdd(&delta[i].y, di.y);
  atomicAdd(&delta[i].z, di.z);
  atomicAdd(&delta[j].x, -dj.x);
  atomicAdd(&delta[j].y, -dj.y);
  atomicAdd(&delta[j].z, -dj.z);

  atomicExch(changed_flag, 1);
}

__global__ void shake_solute_apply_kernel(
    double3* __restrict__ pos,
    double3* __restrict__ delta,
    unsigned num_atoms) {
  const unsigned a = blockIdx.x * blockDim.x + threadIdx.x;
  if (a >= num_atoms) return;
  pos[a] += delta[a];
  delta[a] = double3{0.0, 0.0, 0.0};
}

} // namespace gpu

void gpu::launch_shake_solute_round(
    const double3* pos,
    const double3* old_pos,
    const gpu::ShakeConstraint* constraints,
    unsigned num_constraints,
    const double* inv_mass,
    double tolerance,
    math::boundary_enum boundary,
    math::Box box,
    double dt2,
    double3* delta,
    double3* constraint_force,
    double* virial,
    int* changed_flag,
    int* error_flag,
    cudaStream_t stream) {

  if (num_constraints == 0) return;

  const unsigned threads = 128;
  const unsigned blocks = (num_constraints + threads - 1) / threads;

  switch (boundary) {
    case math::vacuum:
      gpu::shake_solute_round_kernel<math::vacuum><<<blocks, threads, 0, stream>>>(
          pos, old_pos, constraints, num_constraints, inv_mass, tolerance,
          gpu::Periodicity<math::vacuum>(box), dt2, delta, constraint_force,
          virial, changed_flag, error_flag);
      break;
    case math::rectangular:
      gpu::shake_solute_round_kernel<math::rectangular><<<blocks, threads, 0, stream>>>(
          pos, old_pos, constraints, num_constraints, inv_mass, tolerance,
          gpu::Periodicity<math::rectangular>(box), dt2, delta, constraint_force,
          virial, changed_flag, error_flag);
      break;
    case math::triclinic:
      gpu::shake_solute_round_kernel<math::triclinic><<<blocks, threads, 0, stream>>>(
          pos, old_pos, constraints, num_constraints, inv_mass, tolerance,
          gpu::Periodicity<math::triclinic>(box), dt2, delta, constraint_force,
          virial, changed_flag, error_flag);
      break;
    default:
      break;
  }
}

void gpu::launch_shake_solute_apply(
    double3* pos,
    double3* delta,
    unsigned num_atoms,
    cudaStream_t stream) {

  if (num_atoms == 0) return;

  const unsigned threads = 128;
  const unsigned blocks = (num_atoms + threads - 1) / threads;
  gpu::shake_solute_apply_kernel<<<blocks, threads, 0, stream>>>(pos, delta, num_atoms);
}
