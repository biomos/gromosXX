/**
 * @file omp_nonbonded_interaction.h
 * the non bonded interactions:
 * Lennard-Jones and Coulomb interactions.
 * OpenMP code
 */

#ifndef INCLUDED_OMP_NONBONDED_INTERACTION_H
#define INCLUDED_OMP_NONBONDED_INTERACTION_H

#include "nonbonded_interaction.h"

namespace interaction
{
  /**
   * @class OMP_Nonbonded_Interaction
   * calculates the nonbonded interactions using OpenMP
   */
  class OMP_Nonbonded_Interaction : 
    public Nonbonded_Interaction
  {
  public:    
    /**
     * Constructor.
     */
    OMP_Nonbonded_Interaction(Pairlist_Algorithm *pa);
    /**
     * Destructor.
     */
    virtual ~OMP_Nonbonded_Interaction();
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
    /**
     * size the arrays of storage.
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);

  protected:

  };
  
} // interaction

#endif
