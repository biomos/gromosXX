/**
 * @file twin_range_pairlist_cg.h
 * the twin range pairlist (chargegroup - based cutoff).
 */

#ifndef INCLUDED_TWIN_RANGE_PAIRLIST_CG_H
#define INCLUDED_TWIN_RANGE_PAIRLIST_CG_H

namespace interaction
{
  /**
   * @class twin_range_pairlist_cg
   * a pair of simple_pairlists, the
   * distance calculation between atoms is
   * chargegroup based.
   */
  template<typename t_simulation>
  class twin_range_pairlist_cg :
    public twin_range_pairlist<t_simulation>
  {
  public:
    /**
     * update the pairlists (using rcutl, rcutp).
     */
    void update(t_simulation &sim);
  };
} // interaction

#include "twin_range_pairlist_cg.tcc"

#endif
