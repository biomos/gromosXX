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
#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/math/periodicity.h"
#include "math/gmath.h"

#include "shake_kernels.h"

namespace gpu {

template <math::boundary_enum BOUNDARY>
__global__ void shake_solvent_kernel(
    FPL3_TYPE* __restrict__ pos,
    const FPL3_TYPE* __restrict__ old_pos,
    const gpu::ShakeConstraint* __restrict__ constraints,
    unsigned num_constraints,
    const FPL_TYPE* __restrict__ inv_mass_local,
    unsigned num_atoms_per_molecule,
    unsigned first_atom,
    unsigned num_molecules,
    FPL_TYPE tolerance,
    unsigned max_iterations,
    gpu::Periodicity<BOUNDARY> periodicity,
    FPL_TYPE dt2,
    FPH3_TYPE* __restrict__ constraint_force,
    double* __restrict__ virial,
    int* __restrict__ error_flag) {

  const unsigned mol = blockIdx.x * blockDim.x + threadIdx.x;
  if (mol >= num_molecules) return;

  const unsigned base = first_atom + mol * num_atoms_per_molecule;

  // Cache this molecule's atoms in registers/local memory for the whole
  // Jacobi loop instead of re-reading pos[]/old_pos[] from global memory
  // on every constraint check of every iteration. old_pos never changes
  // within apply(), so it's loaded once and never written back; pos and
  // the constraint-force accumulator are written back exactly once,
  // after convergence, instead of per-update.
  FPL3_TYPE local_pos[gpu::MAX_SHAKE_ATOMS_PER_MOLECULE];
  FPL3_TYPE local_old_pos[gpu::MAX_SHAKE_ATOMS_PER_MOLECULE];
  FPL3_TYPE local_cf[gpu::MAX_SHAKE_ATOMS_PER_MOLECULE];
  for (unsigned a = 0; a < num_atoms_per_molecule; ++a) {
    local_pos[a] = pos[base + a];
    local_old_pos[a] = old_pos[base + a];
    local_cf[a] = make_FPL3(FPL_TYPE(0), FPL_TYPE(0), FPL_TYPE(0));
  }

  bool skip_now[gpu::MAX_SHAKE_ATOMS_PER_MOLECULE];
  bool skip_next[gpu::MAX_SHAKE_ATOMS_PER_MOLECULE];
  for (unsigned a = 0; a < num_atoms_per_molecule; ++a) {
    skip_now[a] = false;
    skip_next[a] = true;
  }

  unsigned iterations = 0;
  bool convergence = false;
  // Virial stays high precision (matches the global accumulator it
  // eventually feeds, and lj_crf_tiles.cu's own low-per-thread/
  // high-at-reduction convention) even though the geometry it's built
  // from is FPL_TYPE.
  double v_local[9] = {0,0,0,0,0,0,0,0,0};

  while (!convergence) {
    convergence = true;

    for (unsigned c = 0; c < num_constraints; ++c) {
      const unsigned li = constraints[c].i;
      const unsigned lj = constraints[c].j;

      if (skip_now[li] && skip_now[lj]) continue;
      if (inv_mass_local[li] == FPL_TYPE(0) && inv_mass_local[lj] == FPL_TYPE(0)) continue;

      const FPL3_TYPE r = periodicity.nearest_image(local_pos[li], local_pos[lj]);
      const FPL_TYPE dist2 = dot(r, r);
      const FPL_TYPE r0sq = constraints[c].r0sq;
      const FPL_TYPE diff = r0sq - dist2;

      if (fabs(diff) >= r0sq * tolerance * FPL_TYPE(2)) {
        FPL3_TYPE ref_r = periodicity.nearest_image(local_old_pos[li], local_old_pos[lj]);
        const FPL_TYPE sp = dot(ref_r, r);

        // 1e-12, matching math::epsilon (math/gmath.h) -- that constant
        // isn't device-accessible (a plain non-constexpr host global),
        // so the literal is inlined here instead.
        if (sp < r0sq * FPL_TYPE(1.0e-12)) {
          atomicExch(error_flag, 1);
          return;
        }

        const FPL_TYPE lambda = diff / (sp * FPL_TYPE(2) * (inv_mass_local[li] + inv_mass_local[lj]));

        const FPL3_TYPE cons_force = lambda * ref_r;
        local_cf[li] += cons_force;
        local_cf[lj] -= cons_force;

        const double inv_dt2 = static_cast<double>(lambda) / static_cast<double>(dt2);
        v_local[0] += static_cast<double>(ref_r.x) * static_cast<double>(ref_r.x) * inv_dt2; // (0,0)
        v_local[1] += static_cast<double>(ref_r.x) * static_cast<double>(ref_r.y) * inv_dt2; // (0,1)
        v_local[2] += static_cast<double>(ref_r.x) * static_cast<double>(ref_r.z) * inv_dt2; // (0,2)
        v_local[3] += static_cast<double>(ref_r.y) * static_cast<double>(ref_r.x) * inv_dt2; // (1,0)
        v_local[4] += static_cast<double>(ref_r.y) * static_cast<double>(ref_r.y) * inv_dt2; // (1,1)
        v_local[5] += static_cast<double>(ref_r.y) * static_cast<double>(ref_r.z) * inv_dt2; // (1,2)
        v_local[6] += static_cast<double>(ref_r.z) * static_cast<double>(ref_r.x) * inv_dt2; // (2,0)
        v_local[7] += static_cast<double>(ref_r.z) * static_cast<double>(ref_r.y) * inv_dt2; // (2,1)
        v_local[8] += static_cast<double>(ref_r.z) * static_cast<double>(ref_r.z) * inv_dt2; // (2,2)

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
    // local_cf accumulated in FPL_TYPE (per-molecule register/local
    // memory, throughput-critical) -- widen to FPH only at this final
    // write into the shared accumulator, same convention as every
    // other atomicAdd site in this file. Not an atomicAdd itself: each
    // thread owns a disjoint [base, base+num_atoms_per_molecule) range
    // (one molecule per thread), so plain += is race-free here.
    constraint_force[base + a].x += static_cast<double>(local_cf[a].x);
    constraint_force[base + a].y += static_cast<double>(local_cf[a].y);
    constraint_force[base + a].z += static_cast<double>(local_cf[a].z);
  }

  // Accumulated locally across every constraint/iteration this thread
  // touched, then flushed with one atomicAdd per component instead of
  // one per constraint update -- thousands of molecule-threads otherwise
  // hammer the same 9 global addresses every iteration, serializing the
  // whole kernel.
  for (unsigned k = 0; k < 9; ++k) {
    if (v_local[k] != 0.0) atomicAdd(&virial[k], v_local[k]);
  }
}

} // namespace gpu

void gpu::launch_shake_solvent(
    FPL3_TYPE* pos,
    const FPL3_TYPE* old_pos,
    const gpu::ShakeConstraint* constraints,
    unsigned num_constraints,
    const FPL_TYPE* inv_mass_local,
    unsigned num_atoms_per_molecule,
    unsigned first_atom,
    unsigned num_molecules,
    FPL_TYPE tolerance,
    unsigned max_iterations,
    math::boundary_enum boundary,
    math::Box box,
    FPL_TYPE dt2,
    FPH3_TYPE* constraint_force,
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
    const FPL3_TYPE* __restrict__ pos,
    const FPL3_TYPE* __restrict__ old_pos,
    const gpu::ShakeConstraint* __restrict__ constraints,
    unsigned num_constraints,
    const FPL_TYPE* __restrict__ inv_mass,
    FPL_TYPE tolerance,
    gpu::Periodicity<BOUNDARY> periodicity,
    FPL_TYPE dt2,
    FPL3_TYPE* __restrict__ delta,
    FPH3_TYPE* __restrict__ constraint_force,
    double* __restrict__ virial,
    int* __restrict__ changed_flag,
    int* __restrict__ error_flag) {

  const unsigned c = blockIdx.x * blockDim.x + threadIdx.x;
  if (c >= num_constraints) return;

  const unsigned i = constraints[c].i;
  const unsigned j = constraints[c].j;

  const FPL3_TYPE r = periodicity.nearest_image(pos[i], pos[j]);
  const FPL_TYPE dist2 = dot(r, r);
  const FPL_TYPE r0sq = constraints[c].r0sq;
  const FPL_TYPE diff = r0sq - dist2;

  if (fabs(diff) < r0sq * tolerance * FPL_TYPE(2)) return;

  FPL3_TYPE ref_r = periodicity.nearest_image(old_pos[i], old_pos[j]);
  const FPL_TYPE sp = dot(ref_r, r);

  // 1e-12, matching math::epsilon -- see shake_solvent_kernel's comment.
  if (sp < r0sq * FPL_TYPE(1.0e-12)) {
    atomicExch(error_flag, 1);
    return;
  }

  const FPL_TYPE lambda = diff / (sp * FPL_TYPE(2) * (inv_mass[i] + inv_mass[j]));

  const FPL3_TYPE cons_force = lambda * ref_r;
  atomicAdd(&constraint_force[i].x, static_cast<double>(cons_force.x));
  atomicAdd(&constraint_force[i].y, static_cast<double>(cons_force.y));
  atomicAdd(&constraint_force[i].z, static_cast<double>(cons_force.z));
  atomicAdd(&constraint_force[j].x, -static_cast<double>(cons_force.x));
  atomicAdd(&constraint_force[j].y, -static_cast<double>(cons_force.y));
  atomicAdd(&constraint_force[j].z, -static_cast<double>(cons_force.z));

  const double inv_dt2 = static_cast<double>(lambda) / static_cast<double>(dt2);
  atomicAdd(&virial[0], static_cast<double>(ref_r.x) * static_cast<double>(ref_r.x) * inv_dt2);
  atomicAdd(&virial[1], static_cast<double>(ref_r.x) * static_cast<double>(ref_r.y) * inv_dt2);
  atomicAdd(&virial[2], static_cast<double>(ref_r.x) * static_cast<double>(ref_r.z) * inv_dt2);
  atomicAdd(&virial[3], static_cast<double>(ref_r.y) * static_cast<double>(ref_r.x) * inv_dt2);
  atomicAdd(&virial[4], static_cast<double>(ref_r.y) * static_cast<double>(ref_r.y) * inv_dt2);
  atomicAdd(&virial[5], static_cast<double>(ref_r.y) * static_cast<double>(ref_r.z) * inv_dt2);
  atomicAdd(&virial[6], static_cast<double>(ref_r.z) * static_cast<double>(ref_r.x) * inv_dt2);
  atomicAdd(&virial[7], static_cast<double>(ref_r.z) * static_cast<double>(ref_r.y) * inv_dt2);
  atomicAdd(&virial[8], static_cast<double>(ref_r.z) * static_cast<double>(ref_r.z) * inv_dt2);

  ref_r *= lambda;
  const FPL3_TYPE di = ref_r * inv_mass[i];
  const FPL3_TYPE dj = ref_r * inv_mass[j];
  atomicAdd(&delta[i].x, di.x);
  atomicAdd(&delta[i].y, di.y);
  atomicAdd(&delta[i].z, di.z);
  atomicAdd(&delta[j].x, -dj.x);
  atomicAdd(&delta[j].y, -dj.y);
  atomicAdd(&delta[j].z, -dj.z);

  atomicExch(changed_flag, 1);
}

__global__ void shake_solute_apply_kernel(
    FPL3_TYPE* __restrict__ pos,
    FPL3_TYPE* __restrict__ delta,
    unsigned num_atoms) {
  const unsigned a = blockIdx.x * blockDim.x + threadIdx.x;
  if (a >= num_atoms) return;
  pos[a] += delta[a];
  delta[a] = make_FPL3(FPL_TYPE(0), FPL_TYPE(0), FPL_TYPE(0));
}

} // namespace gpu

void gpu::launch_shake_solute_round(
    const FPL3_TYPE* pos,
    const FPL3_TYPE* old_pos,
    const gpu::ShakeConstraint* constraints,
    unsigned num_constraints,
    const FPL_TYPE* inv_mass,
    FPL_TYPE tolerance,
    math::boundary_enum boundary,
    math::Box box,
    FPL_TYPE dt2,
    FPL3_TYPE* delta,
    FPH3_TYPE* constraint_force,
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
    FPL3_TYPE* pos,
    FPL3_TYPE* delta,
    unsigned num_atoms,
    cudaStream_t stream) {

  if (num_atoms == 0) return;

  const unsigned threads = 128;
  const unsigned blocks = (num_atoms + threads - 1) / threads;
  gpu::shake_solute_apply_kernel<<<blocks, threads, 0, stream>>>(pos, delta, num_atoms);
}
