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
   * The tile holding the interaction pairs dimension 8x4 (based on the mask size)
   */
  struct Interaction_Tile {
      /**
       * Start index of atom/charge group group X (row tile)
       */
      int start_x;
      
      /**
       * Start index of atom/charge group group X (row tile)
       */
      int start_y;

      /**
       * Unique tile ID (do we need it?)
       */
      int index;

      // possibly store other information here instead
      
        // uint16_t flags;      // bit-packed flags: isDiagonal, exclusions, etc.
        // uint8_t atomsX;      // Atom count in X block (optional)
        // uint8_t atomsY;      // Atom count in Y block (optional)

      /**
       * Bitmask of excluded atoms (not interacting, excluded, custom interactions etc.)
       */
      int mask;

      __device__ __host__ Interaction_Tile() : start_x(0), start_y(0), index(0), mask(0) {}

      __device__ __host__ Interaction_Tile(int x, int y, int idx, int mask = 0)
          : start_x(x), start_y(y), index(idx), mask(mask) {}
  };

  /**
   * assert the size is 16 (for efficient alignment)
   */
  static_assert(sizeof(Interaction_Tile) = 16);

  class TileVector {
  public:
      __host__ TileVector(size_t capacity = 0)
          : m_data(nullptr), m_size(nullptr), m_capacity(0), m_overflow(nullptr)
      {
          if (capacity > 0) {
              allocate(capacity);
          }
      }

      __host__ ~TileVector() {
          if (m_data) cudaFree(m_data);
          if (m_size) cudaFree(m_size);
          if (m_overflow) cudaFree(m_overflow);
      }

      __host__ void allocate(size_t capacity) {
          m_capacity = capacity;

          // Allocate Unified Memory
          cudaMallocManaged(&m_data, sizeof(Interaction_Tile) * capacity);
          cudaMallocManaged(&m_size, sizeof(int));
          cudaMallocManaged(&m_overflow, sizeof(bool));

          *_size = 0;
          *_overflow = false;
      }

      __device__ __host__ Interaction_Tile& operator[](size_t i) {
          assert(i < _capacity);
          return _data[i];
      }

      __device__ __host__ const Interaction_Tile& operator[](size_t i) const {
          assert(i < _capacity);
          return _data[i];
      }

      /// Add from device using atomicAdd. Sets overflow flag if exceeded.
      __device__ bool push_back(const Interaction_Tile& tile) {
          int i = atomicAdd(_size, 1);
          if (i >= _capacity) {
              *_overflow = true;  // Signal overflow
              return false;
          }
          _data[i] = tile;
          return true;
      }

      __host__ void resize_host(int new_size) {
          assert(new_size <= _capacity);
          *_size = new_size;
      }

      __device__ __host__ int size() const {
          return *_size;
      }

      __device__ __host__ int capacity() const {
          return _capacity;
      }

      __device__ __host__ Interaction_Tile* data() {
          return _data;
      }

      __device__ __host__ const Interaction_Tile* data() const {
          return _data;
      }

      __host__ bool was_overflown() const {
          return *_overflow;
      }

      __host__ void reset_overflow() {
          if (_overflow) *_overflow = false;
      }

  private:
      Interaction_Tile* _data;
      int* _size;
      int _capacity;
      bool* _overflow;
  };

    /** 
   * @struct TilesContainer
   * holds a set of interacting tiles.
   */
  struct TilesContainer {
    /**
     * resizes all tiles to length 
     */
    inline void resize(unsigned int length) {
      solute_short.resize(length);
      solute_long.resize(length);
      solvent_short.resize(length);
      solvent_long.resize(length);
    }
    
    /**
     * reserve some space 
     */
    inline void reserve(unsigned int pairs) {
      unsigned int n = size();
      for(unsigned int i = 0; i < n; ++i) {
        solute_short[i].reserve(pairs);
        solute_long[i].reserve(pairs);
        solvent_short[i].reserve(pairs);
        solvent_long[i].reserve(pairs);
      }
    }
    
    /** 
     * clears all pairlists
     */
    inline void clear() {
      for(unsigned int i = 0; i < solute_short.size(); ++i) {
        solute_short[i].clear();
        solute_long[i].clear();
        solvent_short[i].clear();
        solvent_long[i].clear();
      }
    }
    /**
     * gives size of pairlists
     */
    inline unsigned int size() const {
      assert(solute_short.size() == solute_long.size() && 
             solute_short.size() == solvent_short.size() &&
             solute_short.size() == solvent_long.size());
      return solute_short.size();
    }
    /**
     * shortrange pairlists that holds: solute-, solute-solute, solute-solvent pairs
     */
    Pairlist solute_short;
    /**
     * longrange pairlists that holds: solute-, solute-solute, solute-solvent pairs
     */
    Pairlist solute_long;
    /**
     * shortrange pairlists that holds: solvent-solvent pairs
     */
    Pairlist solvent_short;
    /**
     * longrange pairlists that holds: solvent-solvent pairs
     */
    Pairlist solvent_long;   
  };

}
