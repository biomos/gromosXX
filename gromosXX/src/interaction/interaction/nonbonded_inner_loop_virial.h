/**
 * @file nonbonded_inner_loop_virial.h
 * inner loop class of the nonbonded routines.
 */

#ifndef INCLUDED_NONBONDED_INNER_LOOP_VIRIAL_H
#define INCLUDED_NONBONDED_INNER_LOOP_VIRIAL_H

namespace interaction
{
  /**
   * @class Nonbonded_Inner_Loop_Virial
   * standard non bonded inner loop.
   * with virial calculation.
   */
  template<typename t_simulation, typename t_storage>
  class Nonbonded_Inner_Loop_Virial
    : Nonbonded_Inner_Loop<t_simulation, t_storage>
  {
  public:
    /**
     * Constructor
     */
    Nonbonded_Inner_Loop_Virial(Nonbonded_Base &base, t_storage &storage);
    

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
    
  };
  
} // interaction

#include "nonbonded_inner_loop_virial.tcc"

#endif

  
    
