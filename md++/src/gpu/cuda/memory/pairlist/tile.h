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
#include "gpu/cuda/cuhostdevice.h"

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
      /**
       * @brief Non-owning, trivially-copyable handle to a TileVecT's
       * device-accessible state -- pass THIS into kernels, never the
       * owning TileVecT itself.
       *
       * TileVecT owns its memory (its destructor calls cudaFree). Kernel
       * launch syntax `kernel<<<...>>>(..., some_tile_vec)` passing a
       * TileVecT BY VALUE constructs a host-side temporary (a shallow
       * copy: same m_data/m_size/m_overflow pointers) whose destructor
       * runs synchronously right after the (asynchronous) launch is
       * enqueued -- freeing the shared buffers while the kernel may still
       * be running or about to run, and while later host code still holds
       * the "original" TileVecT with now-dangling pointers. This was a
       * real, silent bug (found via the pairlist-equivalence test,
       * TILE_PAIRLIST_DESIGN.md §5/§6): reads after such a launch could
       * appear correct for a while (freed memory not yet reused) and then
       * turn to garbage/zero once something else reused the address.
       * View has no destructor, so copying it (including implicitly, as
       * a by-value kernel parameter) is always safe.
       */
      class View {
      public:
          HOSTDEVICE View() : m_data(nullptr), m_size(nullptr), m_capacity(0), m_overflow(nullptr) {}
          HOSTDEVICE View(TileT* data, unsigned* size, unsigned capacity, bool* overflow)
              : m_data(data), m_size(size), m_capacity(capacity), m_overflow(overflow) {}

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

          __device__ __host__ unsigned size() const { return *m_size; }
          __device__ __host__ unsigned capacity() const { return m_capacity; }
          __device__ __host__ TileT* data() { return m_data; }
          __device__ __host__ const TileT* data() const { return m_data; }
          __host__ bool was_overflown() const { return *m_overflow; }

      private:
          TileT *m_data;
          unsigned *m_size;
          unsigned m_capacity;
          bool *m_overflow;
      };

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

      /**
       * @brief Host-only way to populate a TileVecT directly (e.g. a
       * hand-built test tile): both `push_back` overloads are
       * `__device__`-only (atomic-append from a kernel), so there is
       * otherwise no way to set `size()` from host code -- `operator[]`
       * can write tile contents, but `size()` would stay 0 without this.
       * Unlike `push_back`, does not check for overflow: the caller is
       * responsible for `new_size <= capacity()` (use `reserve()` first
       * if needed).
       */
      __host__ void resize(unsigned new_size) {
          reserve(new_size);
          *m_size = new_size;
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

      /**
       * @brief Non-owning handle for kernel parameters -- see View's own
       * doc comment for why this exists and why kernels must take View,
       * never TileVecT itself, by value.
       */
      __host__ View view() const { return View(m_data, m_size, m_capacity, m_overflow); }

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
