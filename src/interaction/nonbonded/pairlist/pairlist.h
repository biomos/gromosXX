/**
 * @file pairlist.h
 * the pairlist class.
 */

#ifndef INCLUDED_PAIRLIST_H
#define INCLUDED_PAIRLIST_H

namespace interaction
{
  /**
   * @class Pairlist
   * holds a Pairlist.
   * very easy implementation that just uses standard vectors.
   */
  class Pairlist :
    public std::vector< std::vector<unsigned int> >
  {
  public:
    /**
     * Constructor.
     */
    Pairlist() {}

  protected:
    
  };
  
  /**
   * sort and print the pairlist.
   */
  std::ostream & 
  operator<<(std::ostream &os, Pairlist &pl);
  
  /** 
   * @struct ParilistContainer
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

#endif
