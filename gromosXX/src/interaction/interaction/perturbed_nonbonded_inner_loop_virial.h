/**
 * @file perturbed_nonbonded_inner_loop_virial.h
 * inner loop class of the perturbed nonbonded routines with virial
 * calculation.
 */

#ifndef INCLUDED_PERTURBED_NONBONDED_INNER_LOOP_VIRIAL_H
#define INCLUDED_PERTURBED_NONBONDED_INNER_LOOP_VIRIAL_H

namespace interaction
{
  /**
   * @class Perturbed_Nonbonded_Inner_Loop_Virial
   * perturbed non bonded inner loop.
   * with virial calculation.
   */
  template<typename t_simulation, typename t_storage>
  class Perturbed_Nonbonded_Inner_Loop_Virial
    : public Nonbonded_Inner_Loop_Virial<t_simulation, t_storage>
  {
  public:
    /**
     * Constructor
     */
    Perturbed_Nonbonded_Inner_Loop_Virial(Nonbonded_Base &base, 
					  t_storage &storage);
    

    /**
     * perturbed interaction
     */
    void perturbed_interaction_inner_loop(t_simulation &sim,
					  size_t const i, size_t const j);
    /**
     * perturbed 1-4 interaction
     */
    void perturbed_one_four_interaction_inner_loop(t_simulation &sim,
					   size_t const i, size_t const j);
    
    
  };
  
} // interaction

#include "perturbed_nonbonded_inner_loop_virial.tcc"

#endif

  
    
