/**
 * @file simple_pairlist.h
 * a simple pairlist class.
 */

#ifndef INCLUDED_SIMPLE_PAIRLIST_H
#define INCLUDED_SIMPLE_PAIRLIST_H

#ifndef INCLUDED_UTILITY
// pair<A,B>
#include <utility>
#define INCLUDED_UTILITY
#endif

namespace interaction
{
  typedef short unsigned int t_pl_index;
  typedef std::vector<t_pl_index> t_pl_row;
  typedef std::vector<t_pl_row> t_pl_matrix;
  
  /**
   * @class simple_pairlist
   * creates a pairlist in the most simple way.
   *
   * @todo add cryptic code somewhere to minimize memory 
   * allocation/deallocation in the update() method.
   */
  template<typename t_simulation>
  class simple_pairlist : 
    public t_pl_matrix
  {
  public:

    class iterator {
    public:
      /**
       * Construct a new iterator over simple_pairlist
       * starting at row(ii), column(jj).
       */
      iterator(t_pl_matrix &pl, 
               int ii = 0,
               int jj = 0);
      /**
       * Go to the next non-empty pairlist entry.
       */
      void operator++();
      /**
       * equality operator.
       */
      bool operator==(iterator it);
      /**
       * unequality operator.
       */
      bool operator!=(iterator it);
      /**
       * Return the value of the pairlist entry.
       */
      t_pl_index& operator*() { return m_pairlist[i()][j()]; }
      /**
       * Return the row index.
       */
      t_pl_index& i() { return m_i; }
      /**
       * Return the column index.
       * @TODO check if properly updated or remove.
       * *it works...
       */
      t_pl_index& j() { return m_j; }
      /**
       * Pairlist accessor.
       */
      t_pl_matrix& pairlist() { return m_pairlist; }
    protected:
      /**
       * reference to the pairlist we operate over.
       */
      t_pl_matrix &m_pairlist;
      /**
       * index of atom i at current position.
       */
      t_pl_index m_i;
      /**
       * index of atom j at current position.
       * @TODO check if properly updated.
       */
      t_pl_index m_j;
    };

    /**
     * Build up a pairlist. Here every particle interacts with 
     * all other particles.
     * To not count double, the pairlist contains only
     * interactions with particles with higher sequence
     * number.
     * Exclusions are checked. But be aware that the exclusions need to
     * be added to the topology (all_exclusions) vector. Currently
     * this is not done for solvent (as exclusions for solvent are
     * normally handled directly).
     */
    void update(t_simulation &sim);
    /**
     * If possible, return an iterator to the first non-empty 
     * pairlist entry. Otherwise return end().
     */
    typename simple_pairlist<t_simulation>::iterator begin();
    typename simple_pairlist<t_simulation>::iterator end();
  };

  template<typename t_simulation>
  std::ostream& operator<<(
			   std::ostream& os, 
			   class simple_pairlist<t_simulation>& pl);
  
} // interaction

// template methods
#include "simple_pairlist.tcc"

#endif
