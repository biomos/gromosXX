/**
 * @file perturbed_nonbonded_outerloop.h
 * the perturbed non bonded outer loop:
 * Lennard-Jones and Coulomb interactions.
 */

#ifndef INCLUDED_PERTURBED_NONBONDED_OUTERLOOP_H
#define INCLUDED_PERTURBED_NONBONDED_OUTERLOOP_H

namespace interaction
{
  /**
   * @class Perturbed_Nonbonded_Outerloop
   * outerloop for the perturbed nonbonded interactions.
   */
  template<typename t_interaction_spec,
	   typename t_perturbation_details>
  class Perturbed_Nonbonded_Outerloop : 
    public Perturbed_Nonbonded_Innerloop<t_interaction_spec, t_perturbation_details>
  {
  public:    

    typedef math::Periodicity<t_interaction_spec::boundary_type> Periodicity_type;

    /**
     * Constructor.
     * @param sim where to store forces and energies
     * (and virial contribution).
     */
    Perturbed_Nonbonded_Outerloop(Nonbonded_Parameter & nbp);

  protected:
    /**
     * calculate the perturbed interactions.
     */
    void perturbed_lj_crf_outerloop(topology::Topology & topo,
				    configuration::Configuration & conf,
				    simulation::Simulation & sim,
				    std::vector<std::vector<unsigned int> > & pl,
				    Storage & storage);

    /**
     * calculate the perturbed 1,4-interactions.
     */
    void perturbed_one_four_outerloop(topology::Topology & topo,
				      configuration::Configuration & conf,
				      simulation::Simulation & sim,
				      Storage & storage);

    /**
     * calculate the perturbed RF contributions for excluded atoms.
     */
    void perturbed_RF_excluded_outerloop(topology::Topology & topo,
					 configuration::Configuration & conf,
					 simulation::Simulation & sim,
					 Storage & storage);
  };
  
} // interaction

// template methods
#include "perturbed_nonbonded_outerloop.cc"

#endif
