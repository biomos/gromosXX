#ifndef INCLUDED_SIMPLE_PAIRLIST_H
#define INCLUDED_SIMPLE_PAIRLIST_H

/**
 * @file simple_pairlist.h
 * a simple pairlist class.
 */

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
      bool operator==(iterator it);
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
       */
      t_pl_index& j() { return m_j; }
      /**
       * Pairlist accessor.
       */
      t_pl_matrix& pairlist() { return m_pairlist; }
    protected:
      t_pl_matrix &m_pairlist;
      t_pl_index m_i;
      t_pl_index m_j;
    };

    /**
     * Build up a pairlist. Here every particle interacts with 
     * all other particles.
     * To not count double, the pairlist contains only
     * interactions with particles with higher sequence
     * number.
     */
    void update(t_simulation &simu);
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
    class simple_pairlist<t_simulation>& pl
  );


  /**
   * @class twin_range_pairlist
   * a pair of simple_pairlists 
   *
   * @todo add cryptic code somewhere to minimize memory 
   * allocation/deallocation in the update() method.
   * 
   * @todo un-hardcode cutoff-distances and nearest-image.
   */
  template<typename t_simulation>
  class twin_range_pairlist :
  public std::pair< simple_pairlist<t_simulation>, simple_pairlist<t_simulation> >
  {

    public:
      /**
       * the short range pairlist
       */
      simple_pairlist<t_simulation>& short_range() { return first; }
      /**
       * the long range pairlist
       */
      simple_pairlist<t_simulation>& long_range() { return second; }
      /**
       * update the pairlists (using rcutl, rcutp).
       */
      void update(t_simulation &simu);
  };
  
} // interaction

// template methods
#include "simple_pairlist.tcc"

#endif
