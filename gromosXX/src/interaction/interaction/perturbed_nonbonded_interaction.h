/**
 * @file perturbed_nonbonded_interaction.h
 * the perturbed non bonded interactions:
 * Lennard-Jones and Coulomb interactions.
 */

#ifndef INCLUDED_PERTURBED_NONBONDED_INTERACTION_H
#define INCLUDED_PERTURBED_NONBONDED_INTERACTION_H

namespace interaction
{
  /**
   * @class Perturbed_Nonbonded_Interaction
   * calculates the nonbonded interactions including perturbed interactions.
   */
  template<typename t_simulation, typename t_pairlist, typename t_innerloop>
  class Perturbed_Nonbonded_Interaction : 
    public Interaction<t_simulation>
  {
  public:    
    /**
     * Constructor.
     * @param sim where to store forces and energies
     * (and virial contribution).
     * @param nonbonded_interaction which interaction functions to use.
     */
    Perturbed_Nonbonded_Interaction(t_simulation &sim, 
    interaction::Nonbonded_Interaction<t_simulation, t_pairlist, t_innerloop>
      & nonbonded_interaction);
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Nonbonded_Interaction();

    /**
     * calculate the interactions.
     */
    virtual void calculate_interactions(t_simulation &sim);

  protected:
    /**
     * calculate the interactions.
     */
    virtual void do_perturbed_interactions(t_simulation &sim,
					   typename t_pairlist::iterator it, 
					   typename t_pairlist::iterator to);

    /**
     * calculate the 1,4-interactions.
     */
    virtual void do_perturbed_14_interactions(t_simulation &sim);

    /**
     * calculate the RF contributions for excluded atoms.
     */
    virtual void do_perturbed_RF_excluded_interactions(t_simulation &sim);

    /**
     * the nonbonded interaction
     */
    interaction::Nonbonded_Interaction<t_simulation, t_pairlist, t_innerloop> 
      & m_nonbonded_interaction;
    
  };
  
} // interaction

// template methods
#include "perturbed_nonbonded_interaction.tcc"

#endif
