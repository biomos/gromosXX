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
     * perturbed interaction
     */
    void perturbed_interaction_inner_loop(t_simulation &sim,
					  size_t const i, size_t const j);
    /**
     * 1-4 interaction
     */
    void one_four_interaction_inner_loop(t_simulation &sim,
					 size_t const i, size_t const j);

    /**
     * perturbed 1-4 interaction
     */
    void perturbed_one_four_interaction_inner_loop(t_simulation &sim,
					   size_t const i, size_t const j);
    

    /**
     * RF interaction (solute).
     */
    void RF_excluded_interaction_inner_loop(t_simulation &sim,
					    size_t const i);

    /**
     * perturbed RF interaction (solute).
     */
    void perturbed_RF_excluded_interaction_inner_loop(t_simulation &sim,
						      std::map<size_t, simulation::Perturbed_Atom>::const_iterator const & mit);

    /**
     * RF solvent interaction.
     */
    void RF_solvent_interaction_inner_loop(t_simulation &sim,
					   simulation::chargegroup_iterator const & cg_it);
    
    /**
     * perturbed pairs.
     */
    void perturbed_pair_interaction_inner_loop(t_simulation &sim,
					       std::vector<simulation::
					       Perturbed_Atompair>::const_iterator const &it);
 
  protected:
    Nonbonded_Base &m_base;
    t_storage &m_storage;
  };
  
} // interaction

#include "nonbonded_inner_loop.tcc"

#endif

  
    
