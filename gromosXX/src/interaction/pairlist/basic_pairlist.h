/**
 * @file basic_pairlist.h
 * the basic pairlist class.
 */

#ifndef INCLUDED_BASIC_PAIRLIST_H
#define INCLUDED_BASIC_PAIRLIST_H

namespace interaction
{
  /**
   * @class basic_pairlist
   * holds a pairlist and provides an iterator.
   */
  template<typename t_simulation, typename t_pairlist_algorithm>
  class Basic_Pairlist :
    public basic_pairlist_type,
    public t_pairlist_algorithm
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
      basic_pairlist_type &m_pairlist;

    };
    
    /**
     * Constructor.
     */
    Basic_Pairlist(interaction::Nonbonded_Base &base);

    /**
     * @returns an iterator over the pairlist.
     */
    Basic_Pairlist::iterator begin();
    /**
     * @returns an iterator to the end of the pairlist.
     */
    Basic_Pairlist::iterator end();
    
    /**
     * @returns an iterator over the perturbed pairlist.
     */
    Basic_Pairlist::iterator perturbed_begin();
    /**
     * @returns an iterator to the end of the perturbed pairlist.
     */
    Basic_Pairlist::iterator perturbed_end();
    
    /**
     * the perturbed pairlist.
     */
    basic_pairlist_type & perturbed();
    /**
     * const accessor to the perturbed pairlist.
     */
    basic_pairlist_type const & perturbed()const;
    
  protected:
    /**
     * the perturbed pairlist.
     * the nonperturbed one is implemented through inheritance...
     */
    basic_pairlist_type m_perturbed_pairlist;
    
  };

  /**
   * print the pairlist.
   */
  template<typename t_simulation, typename t_pairlist_algorithm>
  std::ostream & 
  operator<<(std::ostream &os, Basic_Pairlist<t_simulation, 
	     t_pairlist_algorithm> &pl);
  
} // interaction

#include "basic_pairlist.tcc"

#endif
