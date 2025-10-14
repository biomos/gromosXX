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
 * @file tile.h
 * the interaction tile implementation.
 */

#pragma once

#include <cuda_runtime.h>

namespace gpu
{
//   #define ALIGN32 alignas(32)
  /**
   * The tile holding the interaction pairs dimension 32x32
   */
  template <unsigned ROWS = 32, unsigned COLS = 32>
  struct Interaction_TileT {
      // static constexpr unsigned ROWS = 32;
      // static constexpr unsigned COLS = 32;
      static_assert(ROWS*COLS % (8 /* bits */ * sizeof(unsigned)) == 0);
      static constexpr unsigned MASK_SIZE = ROWS*COLS/(8 * sizeof(unsigned));
      /**
       * Start index of atom/charge group X (row tile)
       */
      // int start_x; // can be back-calculated from index
      
      /**
       * Start index of atom/charge group Y (row tile)
       */
      // int start_y; // can be back-calculated from index

      /**
       * Unique tile ID (do we need it?)
       */
      unsigned index;

      // possibly store other information here instead
      
        // uint16_t flags;      // bit-packed flags: isDiagonal, exclusions, etc.
        // uint8_t atomsX;      // Atom count in X block (optional)
        // uint8_t atomsY;      // Atom count in Y block (optional)

      /**
       * Bitmask of excluded atoms (not interacting, excluded, custom interactions etc.)
       */
      unsigned mask[MASK_SIZE];

      __device__ __host__ Interaction_TileT() : index(0) {
        // #pragma unroll
        for (unsigned i = 0; i < MASK_SIZE; ++i) mask[i] = 0;
      }

      __device__ __host__ Interaction_TileT(int idx, unsigned imask = 0) : index(idx) {
        // #pragma unroll
          for (unsigned i = 0; i < MASK_SIZE; ++i) mask[i] = imask;
      };

      __device__ __host__ Interaction_TileT(int idx, int imask = 0) : index(idx) {
        // #pragma unroll
          for (unsigned i = 0; i < MASK_SIZE; ++i) mask[i] = imask;
      };
  };
  /**
   * assert the size is 8 (for efficient alignment)
   */
  // static_assert(sizeof(Interaction_Tile) = 8);
  using Interaction_Tile = Interaction_TileT<32,32>;

  /**
   * @brief Implementation of a vector of tiles using cuda unified memory
   * 
   * @tparam TileT 
   */
  template <typename TileT>
  class TileVecT {
  public:
      __host__ TileVecT(size_t capacity = 0)
          : m_data(nullptr), m_size(nullptr), m_capacity(0), m_overflow(nullptr)
      {
          // Allocate in Unified Memory
          cudaMallocManaged(&m_size, sizeof(unsigned));
          cudaMallocManaged(&m_overflow, sizeof(bool));
          if (capacity > 0) {
              allocate(capacity);
          }
      }

      __host__ ~TileVecT() {
          if (m_data) cudaFree(m_data);
          if (m_size) cudaFree(m_size);
          if (m_overflow) cudaFree(m_overflow);
      }

      __host__ void allocate(size_t capacity) {
          // Allocate Unified Memory
          cudaMallocManaged(&m_data, sizeof(TileT) * capacity);
          m_capacity = capacity;
          *m_size = 0;
          *m_overflow = false;
      }

      __host__ void deallocate() {
          if (m_data) cudaFree(m_data);
          m_data = nullptr;
          *m_size = 0;
          m_capacity = 0;
          *m_overflow = false;
      }

      __host__ void reserve(unsigned new_capacity) {
          if (new_capacity > m_capacity) {
            deallocate();
            allocate(new_capacity);
          }
      }

      __device__ __host__ TileT& operator[](size_t i) {
          assert(i < m_capacity);
          return m_data[i];
      }

      __device__ __host__ const TileT& operator[](size_t i) const {
          assert(i < m_capacity);
          return m_data[i];
      }

      /// Add from device using atomicAdd. Sets overflow flag if exceeded.
      __device__ bool push_back(const TileT& tile);

      __device__ __host__ unsigned size() const {
          return *m_size;
      }

      __device__ __host__ unsigned capacity() const {
          return m_capacity;
      }

      __device__ __host__ TileT* data() {
          return m_data;
      }

      __device__ __host__ const TileT* data() const {
          return m_data;
      }

      __host__ bool was_overflown() const {
          return *m_overflow;
      }

      __host__ void reset_overflow() {
          if (m_overflow) *m_overflow = false;
      }

      __host__ void clear() {
          cudaMemset(m_data, 0, m_capacity * sizeof(TileT));
          *m_size = 0;
          *m_overflow = false;
      }

  private:
      /**
       * @brief Tiles - data host/device transparent array in unified memory
       * 
       */
      TileT *m_data;
      /**
       * @brief current size of the array writable by device
       * 
       */
      unsigned *m_size;
      /**
       * @brief total array capacity, passed by value to kernels
       * 
       */
      unsigned m_capacity;
      /**
       * @brief overflow flag writable by device, set to True if m_size > m_capacity
       * 
       */
      bool *m_overflow;
  };

  using TileVec = TileVecT<Interaction_Tile>;

  /** 
   * @struct TileContainer
   * holds a set of interacting tiles.
   */
  template <typename TileVecT>
  struct TileContainerT {
    /**
     * reserve some space 
     */
    inline void reserve(unsigned int num_tiles) {
      solute_short.reserve(num_tiles);
      solute_long.reserve(num_tiles);
      solute_candidates.reserve(num_tiles);
      solvent_short.reserve(num_tiles);
      solvent_long.reserve(num_tiles);
      solvent_candidates.reserve(num_tiles);
    }
    
    /** 
     * clears all pairlists
     */
    inline void clear() {
      solute_short.clear();
      solute_long.clear();
      solute_candidates.clear();
      solvent_short.clear();
      solvent_long.clear();
      solvent_candidates.clear();
    }
    
    /**
     * shortrange pairlists that holds: solute-, solute-solute, solute-solvent pairs
     */
    TileVecT solute_short;
    /**
     * longrange pairlists that holds: solute-, solute-solute, solute-solvent pairs
     */
    TileVecT solute_long;
    /**
     * shortrange pairlists that holds: solvent-solvent pairs
     */
    TileVecT solvent_short;
    /**
     * longrange pairlists that holds: solvent-solvent pairs
     */
    TileVecT solvent_long;
    /**
     * shortrange pairlists that holds: solvent-solvent pairs
     */
    TileVecT solute_candidates;
    /**
     * longrange pairlists that holds: solvent-solvent pairs
     */
    TileVecT solvent_candidates;   
  };
  using TileContainer = TileContainerT<TileVec>;

}
