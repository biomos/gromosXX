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
 * @file block_pairlist.h
 * Fixed 32-wide "cluster" block build and block-pair candidate search for
 * the tile-based GPU pairlist. See src/gpu/TILE_PAIRLIST_DESIGN.md §2-3.
 */

#pragma once

#include "gpu/cuda/memory/topology_struct.h"
#include "gpu/cuda/memory/cuvector.h"

namespace gpu {

  constexpr unsigned BLOCK_SIZE = 32;

  /**
   * @brief Pack (row_block, col_block, col_from_b_array) into
   * Interaction_TileT::index. row/col are limited to 15 bits (32767 blocks,
   * i.e. up to ~1M chargegroups per group) -- ample for v1; revisit the
   * encoding if that's ever not enough.
   */
  HOSTDEVICE unsigned pack_block_index(unsigned row_block, unsigned col_block, bool col_from_b) {
    return (row_block & 0x7FFFu) | ((col_block & 0x7FFFu) << 15) | (col_from_b ? (1u << 30) : 0u);
  }

  HOSTDEVICE void unpack_block_index(unsigned index, unsigned& row_block, unsigned& col_block, bool& col_from_b) {
    row_block  = index & 0x7FFFu;
    col_block  = (index >> 15) & 0x7FFFu;
    col_from_b = (index & (1u << 30)) != 0;
  }

  /**
   * @brief Find the chargegroup that owns atom `a`, via binary search over
   * TopologyView::chargegroup (a sorted offset array, size
   * num_chargegroups+1: chargegroup i covers atoms
   * [chargegroup[i], chargegroup[i+1])). Avoids needing a separate
   * per-atom-to-chargegroup array built/uploaded up front.
   */
  HOSTDEVICE unsigned atom_to_chargegroup(const int* chargegroup_offsets, unsigned num_chargegroups, unsigned atom) {
    unsigned lo = 0, hi = num_chargegroups; // hi exclusive
    while (lo + 1 < hi) {
        const unsigned mid = lo + (hi - lo) / 2;
        if (static_cast<unsigned>(chargegroup_offsets[mid]) <= atom) lo = mid;
        else hi = mid;
    }
    return lo;
  }

  /**
   * @brief atom_sort_key[a] = cg_sort_key[atom_to_chargegroup(topo, a)],
   * one thread per atom. Every atom inherits its owning chargegroup's
   * Morton cell key so atoms of the same chargegroup stay spatially
   * adjacent after sorting (TILE_PAIRLIST_DESIGN.md §3 step 1).
   */
  __global__ void atom_sort_key_kernel(
      const int* chargegroup_offsets,
      unsigned num_chargegroups,
      unsigned num_atoms,
      const unsigned* cg_sort_key,
      unsigned* atom_sort_key);

  /**
   * @brief out[i] = offset + i for i in [0, n) -- the identity
   * permutation reorder() seeds before sorting it by Morton key.
   * Hand-written instead of thrust::sequence(thrust::device, ...):
   * on this codebase's largest real-world input (extended_test/
   * ubiquitin, ~22.7k atoms) thrust/cub's internal PtxVersion() device
   * query (cub/util_device.cuh, called from thrust::sequence's
   * "large problem size" dispatch path) fails with cudaErrorInvalidValue
   * from a plain cudaGetDevice() call, immediately poisoning the CUDA
   * context for everything after -- reproduced consistently on this
   * environment's GPU/CUDA-13.3 combination, confirmed via
   * compute-sanitizer and CUDA_LAUNCH_BLOCKING=1 to be the first CUDA
   * failure of the run (not a leaked earlier error), and confirmed
   * absent for small topologies (aladip, 72 atoms) that never reach
   * that thrust dispatch path at all. See KNOWN_ISSUES.md. A one-line
   * grid-stride kernel needs none of cub's dispatch machinery.
   */
  __global__ void sequence_kernel(unsigned* out, unsigned n, unsigned offset);

  /**
   * @brief One compare-exchange stage of an in-place bitonic sort over
   * `keys`/`values` (ascending by key, `values` carried along), `n` a
   * power of two. Standard textbook bitonic-sort-network kernel: thread
   * `i` compares against `i ^ j`, only acting when `i < i^j` (so each
   * pair is only handled once), direction (ascending/descending)
   * determined by whether bit `k` of `i` is set. The host loop below
   * calls this O(log2(n)^2) times with `k` doubling 2..n and `j`
   * halving `k`/2..1 for each `k` -- see bitonic_sort_by_key()'s body.
   *
   * `keys` are 64-bit composite keys (real 32-bit Morton key in the
   * high bits, original position in the low bits, see
   * build_composite_key_kernel in the .cu) rather than the caller's
   * plain 32-bit keys -- makes the compare-exchange network stable
   * (ties broken by original position, matching thrust::sort_by_key's
   * radix sort, which is stable too) without a separate pass. Not just
   * cosmetic: an unstable tie-break here was observed to occasionally
   * (~1-in-15 runs) reclassify a borderline atom pair differently than
   * the CPU reference between otherwise-identical runs.
   */
  __global__ void bitonic_step_kernel(unsigned long long* keys, unsigned* values,
                                       unsigned n, unsigned j, unsigned k);

  /**
   * @brief Ascending sort of `keys_inout[0,n)` by value, with
   * `values_inout[0,n)` permuted identically (stable: ties broken by
   * original position, see bitonic_step_kernel's doc comment) -- a
   * thrust::sort_by_key replacement backed entirely by hand-written
   * kernels (see sequence_kernel's doc comment above for why:
   * thrust/cub's internal device-capability query breaks for any
   * sufficiently large cub dispatch on this codebase's real-world-scale
   * input, on this environment's GPU/CUDA-13.3 combination, and *every*
   * cub-based sort primitive hits the exact same shared utility
   * function, so swapping `thrust::sort_by_key` for a direct
   * `cub::DeviceRadixSort` call would not have helped -- confirmed by
   * moving the failure from thrust::sequence's dispatch to
   * thrust::sort_by_key's radix-sort dispatch, at the exact same line,
   * when only sequence_kernel had been substituted).
   *
   * Bitonic sort needs a power-of-two length; `scratch_keys`/
   * `scratch_values` are resized here to `next_pow2(n)` and reused
   * across calls (the caller owns them so repeated candidate rebuilds
   * don't reallocate every time). Padding slots get composite key
   * `UINT64_MAX` so they sort to the end and are never copied back --
   * safe because `keys_inout` holds Morton cell indices (bounded by
   * grid extent), whose composite form is never genuinely `UINT64_MAX`.
   *
   * If either scratch buffer needs to grow past its current capacity,
   * `stream` is explicitly synchronized first: the underlying
   * std::vector reallocation copies old contents over on the *host*
   * side, which races against any of this same function's own prior
   * call's kernels (queued on `stream`) still reading/writing that old
   * buffer -- caught the hard way via an intermittent one-atom force
   * mismatch in rf_excluded_gpu.t.cc. Only actually costs anything the
   * first time either buffer grows past current capacity (capacity
   * never shrinks), cheap given this only runs at candidate-rebuild
   * cadence.
   *
   * O(n log^2 n) compare-exchanges, i.e. more kernel launches than a
   * radix sort for the same n -- acceptable here since this only runs
   * at candidate-rebuild cadence (TILE_PAIRLIST_DESIGN.md §10 part 2),
   * not every step, and correctness/independence from cub trumps
   * shaving kernel-launch count for this call site.
   */
  void bitonic_sort_by_key(unsigned* keys_inout, unsigned* values_inout, unsigned n,
                            gpu::cuvector<unsigned long long>& scratch_keys,
                            gpu::cuvector<unsigned>& scratch_values,
                            cudaStream_t stream = nullptr);

  /**
   * @brief Compute a bounding sphere (center + radius) for each fixed-size
   * block of `order` entries.
   *
   * Block b covers order[b*BLOCK_SIZE .. min((b+1)*BLOCK_SIZE, count)) --
   * the last block is naturally shorter when count isn't a multiple of
   * BLOCK_SIZE, no padding sentinel needed since the valid range is
   * computed here, not looked up per-entry.
   *
   * @param order permutation: block-sorted position -> chargegroup index
   * @param count number of valid entries in `order`
   * @param num_blocks ceil(count / BLOCK_SIZE)
   * @param cg_cog per-chargegroup representative (box-wrapped) position,
   *   indexed by the *original* chargegroup index (i.e. order[k]): true
   *   cog for solute, first-atom position for solvent (see
   *   Periodicity::prepare_chargegroup).
   */
  __global__ void compute_block_bounds_kernel(
      const unsigned* order,
      unsigned count,
      unsigned num_blocks,
      math::CuVArray::View cg_cog,
      FPL3_TYPE* block_center,
      FPL_TYPE* block_radius);

  /**
   * @brief Test every block-pair's bounding-sphere distance against
   * `cutoff` and push survivors as candidate tiles.
   *
   * O(num_blocks_a * num_blocks_b) per call -- no neighbor-cell pruning of
   * block pairs yet. Deliberate v1 simplification (TILE_PAIRLIST_DESIGN.md):
   * correctness-preserving (a superset test is always safe, it just costs
   * more comparisons), revisit for performance once the pipeline is
   * validated against Standard_Pairlist_Algorithm.
   *
   * Pass order_a == order_b (and matching bounds arrays) for the
   * self-pair case (solute-solute, solvent-solvent): only (bi <= bj) is
   * tested then, to avoid pushing both (i,j) and (j,i). For the
   * solute-solvent case, pass the solute arrays as "a" and solvent as "b";
   * every (bi, bj) pair is tested and the tile's col_from_b flag records
   * that column indices are solvent-order indices, not solute-order.
   */
  template <math::boundary_enum BOUNDARY>
  __global__ void find_block_candidates_kernel(
      unsigned num_blocks_a,
      unsigned num_blocks_b,
      const FPL3_TYPE* block_center_a, const FPL_TYPE* block_radius_a,
      const FPL3_TYPE* block_center_b, const FPL_TYPE* block_radius_b,
      bool self_pairs,
      Periodicity<BOUNDARY> periodicity,
      FPL_TYPE cutoff,
      TileVecT<Interaction_Tile>::View candidates);

  /**
   * @brief Exclusion + short/long classification (TILE_PAIRLIST_DESIGN.md
   * §3 step 5). One CUDA block per candidate tile; blockDim = (32,32),
   * one warp per tile row so __ballot_sync builds each row's mask word
   * without atomics.
   *
   * `row_order`/`row_count` are always the tile's row-side atom order
   * (and, when a tile's col_from_b flag is false, also its column side --
   * see find_block_candidates_kernel's self_pairs convention).
   * `col_other_order`/`col_other_count` are only used when col_from_b is
   * true (the solute-solvent case): pass the solvent atom order/count
   * when processing solute_candidates, and {nullptr, 0} when processing
   * solvent_candidates (col_from_b is never set on those tiles).
   *
   * ATOMIC_CUTOFF (TILE_PAIRLIST_DESIGN.md §4.2/step 7) selects the
   * distance-test input, same kernel/tiles otherwise: false uses each
   * pair's owning chargegroups' cog (`cg_cog`, via `topo.chargegroup` for
   * the atom->chargegroup lookup); true uses the atoms' own positions
   * (`pos`) directly, matching Standard_Pairlist_Algorithm::update_atomic.
   * Same-chargegroup pairs are handled differently by mode too, matching
   * their respective CPU references exactly:
   *  - solvent same-chargegroup (same molecule) pairs are always skipped
   *    in both modes -- structural (loop-bounds) in the CPU atomic code,
   *    not exclusion-list-based, since solvent atoms have no CSR entries
   *    (see TopologyView::excl_ptr's doc comment).
   *  - solute same-chargegroup pairs: chargegroup-cutoff mode assumes
   *    they're always in range (no distance test, matches
   *    Standard_Pairlist_Algorithm::_update_cg's direct push to
   *    solute_short after only an exclusion check); atomic-cutoff mode
   *    does a real atom-atom distance test like any other pair (matches
   *    update_atomic, which never special-cases same-chargegroup pairs).
   */
  template <bool ATOMIC_CUTOFF, math::boundary_enum BOUNDARY>
  __global__ void classify_tiles_kernel(
      TileVecT<Interaction_Tile>::View candidates,
      const unsigned* row_order, unsigned row_count,
      const unsigned* col_other_order, unsigned col_other_count,
      math::CuVArray::View cg_cog,
      math::CuVArray::View pos,
      Topology::View topo,
      Periodicity<BOUNDARY> periodicity,
      FPL_TYPE cutoff_short2, FPL_TYPE cutoff_long2,
      TileVecT<Interaction_Tile>::View out_short,
      TileVecT<Interaction_Tile>::View out_long);

}
