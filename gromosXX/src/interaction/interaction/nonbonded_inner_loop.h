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
  template<typename t_simulation, typename t_storage>
  class Nonbonded_Inner_Loop
  {
  public:
    /**
     * Constructor.
     */
    Nonbonded_Inner_Loop(Nonbonded_Base &base, t_storage &storage);
    
    /**
     * interaction
     */
    void interaction_inner_loop(t_simulation const &sim,
				size_t const i, size_t const j);

    /**
     * 1-4 interaction
     */
    void one_four_interaction_inner_loop(t_simulation &sim,
					 size_t const i, size_t const j);
    
  protected:
    Nonbonded_Base &m_base;
    t_storage &m_storage;
  };
  
} // interaction

#include "nonbonded_inner_loop.tcc"

#endif

  
    
