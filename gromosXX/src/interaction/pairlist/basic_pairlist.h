/**
 * @file basic_pairlist.h
 * the basic pairlist class.
 */

#ifndef INCLUDED_BASIC_PAIRLIST_H
#define INCLUDED_BASIC_PAIRLIST_H

namespace interaction
{
  typedef std::vector<std::vector<unsigned int> > basic_pairlist_type;
  
  /**
   * @class basic_pairlist
   * holds a pairlist and provides an iterator.
   */
  template<typename t_simulation>
  class Basic_Pairlist :
    public basic_pairlist_type
  {
  public:
    class iterator
    {
    public:
      /**
       * Constructor.
       */
      iterator(basic_pairlist_type &pl);
      
      /**
       * next entry.
       */
      void operator++();
      /**
       * equality operator.
       * only compares i! (for speed)
       */
      bool operator==(iterator &it);
      /**
       * unequality operator.
       * only compares i! (for speed)
       */
      bool operator!=(iterator &it);
      /**
       * the pair: i
       */
      unsigned int i();
      /**
       * the pair: j
       */
      unsigned int j();
      
      void row(unsigned int i);
      
    protected:
      std::vector<std::vector<unsigned int> >::iterator m_i;
      std::vector<unsigned int>::iterator m_j;
      basic_pairlist_type &m_pairlist;
    };
    
    Basic_Pairlist::iterator begin();
    Basic_Pairlist::iterator end();
  };

  template<typename t_simulation>
  std::ostream & 
  operator<<(std::ostream &os, Basic_Pairlist<t_simulation> &pl);
  
} // interaction

#include "basic_pairlist.tcc"

#endif
