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
   * calculates the perturbed nonbonded interactions.
   */
  template<typename t_simulation, typename t_nonbonded_spec>
  class Perturbed_Nonbonded_Interaction : 
    public Nonbonded_Interaction<t_simulation, t_nonbonded_spec>,
    public t_nonbonded_spec::perturbation_filter_type,
    public t_nonbonded_spec::perturbed_nonbonded_innerloop_type
  {
  public:    
    /**
     * Constructor.
     * @param sim where to store forces and energies
     * (and virial contribution).
     */
    Perturbed_Nonbonded_Interaction(t_simulation &sim);
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Nonbonded_Interaction();

    /**
     * calculate the interactions.
     */
    virtual void calculate_interactions(t_simulation &sim);

    /**
     * add a shortrange interaction.
     */
    void add_shortrange_pair(t_simulation const &sim,
				     size_t const i, size_t const j);
    /**
     * add a longrange interaction.
     */
    void add_longrange_pair(t_simulation & sim,
				    size_t const i, size_t const j);

  protected:
    /**
     * calculate the perturbed interactions.
     */
    virtual void do_perturbed_interactions(t_simulation &sim,
					   Pairlist::iterator it, 
					   Pairlist::iterator to);
    /**
     * calculate the perturbed 1,4-interactions.
     */
    virtual void do_perturbed_14_interactions(t_simulation &sim);

    /**
     * calculate the perturbed RF contributions for excluded atoms.
     */
    virtual void do_perturbed_RF_excluded_interactions(t_simulation &sim);

    /**
     * calculate the perturbed pair contributions.
     */
    void do_perturbed_pair_interactions(t_simulation &sim);
  };
  
} // interaction

// template methods
#include "perturbed_nonbonded_interaction.tcc"

#endif
