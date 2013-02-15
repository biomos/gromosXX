/**
 * @file pairlist.h
 * the pairlist class.
 */

#ifndef INCLUDED_PAIRLIST_H
#define INCLUDED_PAIRLIST_H

#ifdef XXCUDA
#include "../../../cuda/gromos_cuda.h"
#endif

namespace interaction
{

#ifdef XXCUDA
  /**
   * @class Pairlist
   * holds a Pairlist
   * implementation that uses gcuda::TwoDArray
   */
  class Pairlist
  {
  public:

      unsigned int width;
      unsigned int height;
      typedef gcuda::TwoDArray_iterator<unsigned int>  iterator;

    /**
     * Constructor.
     */
    Pairlist(unsigned int i, unsigned int j) : vector(i, j) {

        width = i;
        height = j;

    }
     /**
     * Deconstructor.
     */
    ~Pairlist() {}

      /*
      *
      * Methods
      *
      */
   void clear() {
      for (unsigned int i = 0; i < vector.dim_i; i++) {
                vector.elements[i] = 0;
            }
   }

    unsigned int size() {
        return vector.dim_i;
    }

    unsigned int size(unsigned int i) {
        return vector.elements[i];
    }

    void push_back(unsigned int position, unsigned int value) {
         vector.ptr[position][vector.elements[position]] = value;
            vector.elements[position] = vector.elements[position] + 1;
    }

    Pairlist::iterator begin(unsigned int i) {
        return iterator(&vector.ptr[i][0]);

    }

    Pairlist::iterator end(unsigned int i) {
       return  iterator(&vector.ptr[i][vector.elements[i]]);
    }

    gcuda::TwoDArray<unsigned int>& container() {
         return vector;
     }

  protected:
     gcuda::TwoDArray<unsigned int> vector;

  };

#else
   /**
   * @class Pairlist
   * holds a Pairlist.
   * very easy implementation that just uses standard vectors.
   */
  class Pairlist
  {
  public:

      unsigned int width;
      unsigned int height;
      typedef std::vector<unsigned int>::iterator iterator;

    /**
     * Constructor.
     */
    Pairlist(unsigned int i, unsigned int j) : vector(i, std::vector<unsigned int>(j,0)) {
        vector.resize(i);
        unsigned int n = vector.size();
        for(unsigned int k = 0; k < n; ++k) {
            vector[k].reserve(j);
        }

        width = i;
        height = j;
    }
     /**
     * Deconstructor.
     */
    ~Pairlist() {}

      /*
      *
      * Methods
      *
      */
    void clear() {
        for(unsigned int i = 0; i < vector.size(); ++i) {
            vector[i].clear();
        }
    }

    unsigned int size() {
        return vector.size();
    }
    unsigned int size(unsigned int i) {
        return vector[i].size();
    }

    void push_back(unsigned int i, unsigned int j) {
        vector[i].push_back(j);
    }

    Pairlist::iterator begin(unsigned int i) {
        return vector[i].begin();
    }

    Pairlist::iterator end(unsigned int i) {
        return vector[i].end();
    }

     std::vector<std::vector<unsigned int> >& container() {
         return vector;
     }

  protected:
     std::vector<std::vector<unsigned int> > vector;

  };
  #endif
  /**
   * sort and print the pairlist.
   */
  //std::ostream &
  //operator<<(std::ostream &os, Pairlist &pl);
  
    /**
   * @class ParilistContainer
   * holds a set of pairlists.
   */
  // this needs to be a class to be able to be initialize
  //the width and height of the pairlists

  class PairlistContainer {

  public:

      unsigned int dim_i;
      unsigned int dim_j;

    /**
     * Constructor.
     */
     PairlistContainer(unsigned int i, unsigned int j) : solute_short(i, j), solute_long(i, j), solvent_short(i, j), solvent_long(i, j) {

         dim_i = i;
         dim_j = j;

     }
    /**
     * Deconstructor.
     */
     ~PairlistContainer() {}

    /**
     * resizes all pairlists to length
     */
   //deprecated

    /**
     * reserve some space
     */
     // deprecated

    /**
     * clears all pairlists
     */
   void clear() {
        // clears the size matrix, no need to 0 the array!?
        solute_short.clear();
        solute_long.clear();
        solvent_short.clear();
        solvent_long.clear();
    }
    /**
     * gives size of pairlists
     */
    unsigned int size() const {
        return dim_i;
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
