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
      
      void row(size_t i);
      
    protected:
      std::vector<std::vector<size_t> >::iterator m_i;
      std::vector<size_t>::iterator m_j;
      basic_pairlist_type &m_pairlist;

    };
    
    /**
     * Constructor.
     */
    Basic_Pairlist(interaction::Nonbonded_Base &base);

    Basic_Pairlist::iterator begin();
    Basic_Pairlist::iterator end();

    Basic_Pairlist::iterator perturbed_begin();
    Basic_Pairlist::iterator perturbed_end();
    
    basic_pairlist_type & perturbed();
    basic_pairlist_type const & perturbed()const;
    
  protected:
    basic_pairlist_type m_perturbed_pairlist;
    
  };

  template<typename t_simulation, typename t_pairlist_algorithm>
  std::ostream & 
  operator<<(std::ostream &os, Basic_Pairlist<t_simulation, 
	     t_pairlist_algorithm> &pl);
  
} // interaction

#include "basic_pairlist.tcc"

#endif
