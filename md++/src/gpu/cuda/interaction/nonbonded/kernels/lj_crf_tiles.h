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
 * @file lj_crf_tiles.h
 * Host-callable entry point for the LJ + reaction-field force/energy
 * kernel that consumes the tile-based GPU pairlist
 * (`gpu::TileVecT<Interaction_Tile>`), replacing the discarded flat-pair-
 * array kernel (PLAN.md §7/§10 step 8). See
 * src/gpu/TILE_PAIRLIST_DESIGN.md.
 *
 * Deliberately does *not* include `gpu/cuda/math/periodicity.h` or
 * declare the `__global__` kernel template itself: `gpu::Periodicity`
 * isn't host-compilable (its `set_cell_size`/`get_cell` bodies use bare
 * `min`/`max`, which only resolve under `nvcc`), and a `__global__`
 * template declaration is only ever meaningful to `nvcc` anyway. This
 * header is safe to include from plain (non-`nvcc`) `.cc` translation
 * units -- e.g. a correctness test that only needs to launch the kernel,
 * not name its template parameters. `lj_crf_tiles.cu` defines the actual
 * kernel privately and instantiates it directly from
 * `launch_lj_crf_tiles`'s switch below.
 */

#pragma once

#include "gpu/cuda/memory/precision.h"
#include "gpu/cuda/memory/pairlist/tile.h"
#include "gpu/cuda/memory/cuvector.h"
#include "gpu/cuda/interaction/nonbonded/cuda_lj_params.h"
#include "gpu/cuda/interaction/nonbonded/cuda_nb_sim_params.h"
#include "math/gmath.h"
#include "math/box.h"

namespace gpu {

  /**
   * @brief LJ + reaction-field forces and energies for one bucket of
   * classified tiles (e.g. `m_tiles.solute_short`), dispatched to the
   * right `gpu::Periodicity<BOUNDARY>` instantiation for a runtime
   * `math::boundary_enum`.
   *
   * One CUDA block per tile, blockDim = (32, 32) -- one warp per tile row,
   * same layout as `classify_tiles_kernel` (block_pairlist.h), whose
   * `row_order`/`col_other_order` convention this reuses exactly: pass the
   * tile's row-side atom order as `row_order`, and (only consulted when a
   * tile's `col_from_b` flag is set, i.e. a solute-solvent tile)
   * `col_other_order`/`col_other_count` for the column side.
   *
   * Per active (mask bit set) pair, computes exactly
   * `interaction::Nonbonded_Term::lj_crf_interaction`'s default case
   * (`nonbonded_term.cc`: `eps = 0`, `coulomb_scaling = 1` -- no
   * coarse-graining, no 1-4 scaling; those are out of scope for this
   * kernel, same as the rest of the tile pairlist so far):
   *
   *   r        = nearest_image(pos(a1), pos(a2))     (vector a2 -> a1)
   *   e_lj     = (c12/r^12 - c6/r^6)
   *   e_crf    = q * four_pi_eps_i * (1/r - crf_2cut3i*r^2 - crf_cut)
   *   f_scalar = 12*c12/r^14 - 6*c6/r^8
   *              + q * four_pi_eps_i * (1/r^3 + 2*crf_2cut3i)
   *   F_a1 += f_scalar * r,  F_a2 -= f_scalar * r
   *
   * (`crf_cut3i`, used in the CPU force formula, is always `2 *
   * crf_2cut3i` by construction -- see `Nonbonded_Term::init`,
   * `nonbonded_term.cc` -- so `NbSimParams` only needs to carry
   * `crf_2cut3i`, not a separate field.)
   *
   * Forces are accumulated into `force` via `atomicAdd` (one thread per
   * pair, no reduction needed). Energies are bucketed per energy-group
   * pair (`atom_energy_group[a1] * nb.num_energy_groups +
   * atom_energy_group[a2]`) via dynamic-shared-memory atomics (zeroed at
   * the start of the tile, one `atomicAdd` per active pair into the
   * tile's shared bucket, then one `atomicAdd` per bucket -- not per pair
   * -- into `e_lj_total`/`e_crf_total` at the end), accumulated as
   * `double` regardless of `FPL_TYPE` so summing many small
   * contributions doesn't lose precision the way thousands of individual
   * float atomicAdds would. `e_lj_total`/`e_crf_total` must each be sized
   * `nb.num_energy_groups * nb.num_energy_groups`, flattened row-major
   * (matching `configuration::Energy::lj_energy`/`crf_energy`'s
   * `[gi][gj]` shape) -- for the common single-energy-group case this
   * degenerates to exactly one bucket, i.e. one `atomicAdd` per tile,
   * same cost as before.
   *
   * `virial_total` accumulates the atomic virial the same way (shared-
   * memory bucketed, one flush per tile), always exactly 9 elements
   * (flattened row-major 3x3, `b*3+a`) regardless of energy-group count
   * -- the CPU inner loop never buckets virial by energy group either.
   * `virial_total[b*3+a] += r(b) * force(a)` per active pair, the exact
   * `nonbonded_innerloop.cc` formula (`r` = `rvec` below, `force(a)` =
   * `fr(a)`). Always accumulated, independent of whether a virial was
   * actually requested -- matches CPU convention; the caller decides
   * whether to use it.
   *
   * Does not synchronize -- call `cudaDeviceSynchronize()` (or check a
   * stream/event) before reading
   * `force`/`e_lj_total`/`e_crf_total`/`virial_total`, same convention as
   * the pairlist kernels' launch sites (`cuda_pairlist_algorithm_impl.cu`).
   */
  void launch_lj_crf_tiles(
      TileVecT<Interaction_Tile>::View tiles,
      const unsigned* row_order, unsigned row_count,
      const unsigned* col_other_order, unsigned col_other_count,
      math::CuVArray::View pos,
      const int* iac,
      const FPL_TYPE* charge,
      const unsigned* atom_energy_group,
      LJParamView lj,
      NbSimParams nb,
      math::boundary_enum boundary,
      math::Box box,
      FPH3_TYPE* force,
      double* e_lj_total,
      double* e_crf_total,
      double* virial_total,
      cudaStream_t stream = nullptr);

}
