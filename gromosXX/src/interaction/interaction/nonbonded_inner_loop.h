/**
 * @file nonbonded_inner_loop.h
 * inner loop class of the nonbonded routines.
 */

#ifndef INCLUDED_NONBONDED_INNER_LOOP_H
#define INCLUDED_NONBONDED_INNER_LOOP_H

namespace interaction
{
  /**
   * @class Nonbonded_Inner_Loop
   * standard non bonded inner loop.
   * no virial calculation.
   */
  template<typename t_simulation, typename t_pairlist, typename t_storage>
  class Nonbonded_Inner_Loop
  {
  public:
    /**
     * Constructor.
     */
    Nonbonded_Inner_Loop(t_storage &store, Nonbonded_Base &base);
    
    /**
     * interaction
     */
    void do_interaction(t_simulation &sim, typename t_pairlist::iterator &it);

    /**
     * 1-4 interaction
     */
    void do_one_four_interaction(t_simulation &sim, typename t_pairlist::iterator &it);
    
  protected:
    Nonbonded_Base &m_base;
    t_storage &m_storage;
  };
  
} // interaction

#include "nonbonded_inner_loop.tcc"

#endif

  
    
