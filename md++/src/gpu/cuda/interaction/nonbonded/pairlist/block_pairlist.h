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
      TileVecT<Interaction_Tile> candidates);

}
