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
  template<typename t_simulation, 
	   typename t_filter = interaction::Simple_Filter<t_simulation> >
  class twin_range_pairlist :
    public std::pair< simple_pairlist<t_simulation, t_filter>, 
		      simple_pairlist<t_simulation, t_filter> >
  {
  public:
    /**
     * access to the iterator type.
     */
    typedef typename simple_pairlist<t_simulation, t_filter>::iterator iterator;
    
    /**
     * the short range pairlist
     */
    simple_pairlist<t_simulation, t_filter>& short_range() { 
      return std::pair< simple_pairlist<t_simulation, t_filter>, 
	simple_pairlist<t_simulation, t_filter> >::first; }
    /**
     * the long range pairlist
     */
    simple_pairlist<t_simulation, t_filter>& long_range() { 
      return std::pair< simple_pairlist<t_simulation, t_filter>, 
	simple_pairlist<t_simulation, t_filter> >::second; }
    /**
     * update the pairlists (using rcutl, rcutp).
     */
    void update(t_simulation &sim);

  protected:

  };  
  
} // interaction

#include "twin_range_pairlist.tcc"

#endif
