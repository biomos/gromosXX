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
    public std::vector< std::vector<unsigned int> >
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
      iterator(std::vector< std::vector<unsigned int> > &pl);      
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
      /**
       * also the pair: j
       */
      unsigned int operator*();
      /**
       * the row.
       */
      void row(unsigned int i);
      
    protected:
      /** 
       * iterator over the atoms i of atom pairs i - j.
       */
      std::vector<std::vector<unsigned int> >::iterator m_i;
      /**
       * iterator over the atoms j of atom pairs i - j.
       */
      std::vector<unsigned int>::iterator m_j;
      /**
       * reference to the pairlist.
       */
      std::vector<std::vector<unsigned int> > &m_pairlist;

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
  std::ostream & 
  operator<<(std::ostream &os, Pairlist &pl);
  
} // interaction

#include "pairlist.cc"

#endif
