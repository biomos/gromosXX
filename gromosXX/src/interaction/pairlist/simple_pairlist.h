/**
 * @file simple_pairlist.h
 * a simple pairlist class.
 */

#ifndef INCLUDED_SIMPLE_PAIRLIST_H
#define INCLUDED_SIMPLE_PAIRLIST_H

namespace interaction
{
  /**
   * @class simple_pairlist
   * creates a pairlist in the most simple way.
   */
  template<typename t_simulation>
  class simple_pairlist
  {
  public:
    /**
     * @class iterator
     * iterator over the pairlist.
     */
    class iterator
    {
    public:
      iterator(std::vector< std::vector<int> > &pl);
      bool eol();
      void operator++();
      int i();
      int j();
    protected:
      std::vector< std::vector<int> > &m_pairlist;
      int m_atom;
      std::vector< std::vector<int> >::const_iterator m_i;
      std::vector<int>::const_iterator m_j;
    };

    void make_pairlist(t_simulation &simu);
    void print_pairlist(std::ostream &os);
    iterator begin();
  protected:
    void clear_pairlist();
    std::vector<std::vector<int> > m_pairlist;

  };
  
} // interaction

// template methods
#include "simple_pairlist.tcc"

#endif
