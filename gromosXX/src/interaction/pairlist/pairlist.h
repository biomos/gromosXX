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
   * holds a Pairlist and provides an iterator.
   */
  class Pairlist :
    public std::vector< std::vector<size_t> >
  {
  public:
    /**
     * @class iterator
     * provide an iterator over the pairlist.
     */
    class iterator
    {
    public:
      /**
       * Constructor.
       */
      iterator(std::vector< std::vector<size_t> > &pl);      
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
      size_t i();
      /**
       * the pair: j
       */
      size_t j();
      /**
       * also the pair: j
       */
      size_t operator*();
      /**
       * the row.
       */
      void row(size_t i);
      
    protected:
      /** 
       * iterator over the atoms i of atom pairs i - j.
       */
      std::vector<std::vector<size_t> >::iterator m_i;
      /**
       * iterator over the atoms j of atom pairs i - j.
       */
      std::vector<size_t>::iterator m_j;
      /**
       * reference to the pairlist.
       */
      std::vector<std::vector<size_t> > &m_pairlist;

    };
    
    /**
     * Constructor.
     */
    Pairlist();

    /**
     * @returns an iterator over the pairlist.
     */
    Pairlist::iterator begin();
    /**
     * @returns an iterator to the end of the pairlist.
     */
    Pairlist::iterator end();
    
  protected:
    
  };

  /**
   * print the pairlist.
   */
  template<typename t_simulation, typename t_pairlist_algorithm>
  std::ostream & 
  operator<<(std::ostream &os, Pairlist &pl);
  
} // interaction

#include "pairlist.tcc"

#endif
