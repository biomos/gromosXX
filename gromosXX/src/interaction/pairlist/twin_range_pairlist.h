/**
 * @file twin_range_pairlist.h
 * the twin range pairlist (atom - based cutoff).
 */

#ifndef INCLUDED_TWIN_RANGE_PAIRLIST_H
#define INCLUDED_TWIN_RANGE_PAIRLIST_H

namespace interaction
{
  /**
   * @class twin_range_pairlist
   * a pair of simple_pairlists 
   *
   * @todo add cryptic code somewhere to minimize memory 
   * allocation/deallocation in the update() method.
   */
  template<typename t_simulation>
  class twin_range_pairlist :
    public std::pair< simple_pairlist<t_simulation>, simple_pairlist<t_simulation> >
  {
  public:
    /**
     * access to the iterator type.
     */
    typedef typename simple_pairlist<t_simulation>::iterator iterator;
    
    /**
     * the short range pairlist
     */
    simple_pairlist<t_simulation>& short_range() { 
      return std::pair< simple_pairlist<t_simulation>, 
	simple_pairlist<t_simulation> >::first; }
    /**
     * the long range pairlist
     */
    simple_pairlist<t_simulation>& long_range() { 
      return std::pair< simple_pairlist<t_simulation>, 
	simple_pairlist<t_simulation> >::second; }
    /**
     * update the pairlists (using rcutl, rcutp).
     */
    void update(t_simulation &sim);

  };  
  
} // interaction

#include "twin_range_pairlist.tcc"

#endif
