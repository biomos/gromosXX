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
 * @file pairlist.h
 * the CUDA pairlist class.
 */

#pragma once

namespace gpu
{
  /**
   * @class CUDA_Pairlist
   * holds a Pairlist on device memory.
   * very easy implementation that just uses standard vectors.
   */
  class CUDA_Pairlist {
  public:
    /**
     * Constructor.
     */
    CUDA_Pairlist() {}

  protected:
    
  };
  
  /**
   * sort and print the pairlist.
   */
  std::ostream & 
  operator<<(std::ostream &os, Pairlist &pl);
  
  /** 
   * @struct PairlistContainer
   * holds a set of pairlists.
   */
  struct PairlistContainer {
    /**
     * resizes all pairlists to length 
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
  
} // interaction

