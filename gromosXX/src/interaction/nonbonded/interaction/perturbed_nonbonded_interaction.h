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
  template<typename t_interaction_spec>
  class Perturbed_Nonbonded_Interaction : 
    public Nonbonded_Interaction<t_interaction_spec>,
    public t_interaction_spec::perturbation_filter_type,
    public t_interaction_spec::perturbed_nonbonded_innerloop_type
  {
  public:    
    /**
     * Constructor.
     * @param sim where to store forces and energies
     * (and virial contribution).
     */
    Perturbed_Nonbonded_Interaction();
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Nonbonded_Interaction();

    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);

    /**
     * add a shortrange interaction.
     */
    void add_shortrange_pair(topology::Topology & topo,
			     configuration::Configuration & conf,
			     simulation::Simulation & sim,
			     size_t const i, size_t const j);

    /**
     * add a shortrange interaction.
     */
    void add_shortrange_pair(topology::Topology & topo,
			     configuration::Configuration & conf,
			     simulation::Simulation & sim,
			     size_t const i, size_t const j,
			     int pc);

    /**
     * add a longrange interaction.
     */
    void add_longrange_pair(topology::Topology & topo,
			    configuration::Configuration & conf,
			    simulation::Simulation & sim,
			    size_t const i, size_t const j,
			    math::Periodicity<t_interaction_spec::boundary_type>
			    const & periodicity);

    /**
     * add a longrange interaction.
     */
    void add_longrange_pair(topology::Topology & topo,
			    configuration::Configuration & conf,
			    simulation::Simulation & sim,
			    size_t const i, size_t const j,
			    math::Periodicity<t_interaction_spec::boundary_type>
			    const & periodicity, int pc);

  protected:
    /**
     * calculate the perturbed interactions.
     */
    virtual void do_perturbed_interactions(topology::Topology & topo,
					   configuration::Configuration & conf,
					   simulation::Simulation & sim,
					   Pairlist::iterator it, 
					   Pairlist::iterator to);
    /**
     * calculate the perturbed 1,4-interactions.
     */
    virtual void 
    do_perturbed_14_interactions(topology::Topology & topo,
				 configuration::Configuration & conf,
				 simulation::Simulation & sim);
    

    /**
     * calculate the perturbed RF contributions for excluded atoms.
     */
    virtual void 
    do_perturbed_RF_excluded_interactions(topology::Topology & topo,
					  configuration::Configuration & conf,
					  simulation::Simulation & sim);

    /**
     * calculate the perturbed pair contributions.
     */
    void do_perturbed_pair_interactions(topology::Topology & topo,
					configuration::Configuration & conf,
					simulation::Simulation & sim);
    
  };
  
} // interaction

// template methods
#include "perturbed_nonbonded_interaction.tcc"

#endif
