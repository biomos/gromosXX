/**
 * @file nonbonded_set.h
 * the non bonded interactions for a set of atoms:
 * Lennard-Jones and Coulomb interactions
 * perturbed and non-perturbed.
 */

#ifndef INCLUDED_NONBONDED_SET_H
#define INCLUDED_NONBONDED_SET_H

#include "pairlist.h"
#include "storage.h"
#include "nonbonded_outerloop.h"
#include "nonbonded_set_interface.h"

namespace interaction
{
  /**
   * @class Nonbonded_Set
   * calculates the nonbonded interactions.
   */
  class Nonbonded_Set : public Nonbonded_Set_Interface
  {
  public:    
    /**
     * Constructor.
     */
    Nonbonded_Set(Pairlist_Algorithm & pairlist_alg, Nonbonded_Parameter & param,
		  int rank, int num_threads);
    
    /**
     * Destructor
     */
    virtual ~Nonbonded_Set() {}
    
    /**
     * initialize some things
     */
    virtual int init(topology::Topology const & topo,
		     configuration::Configuration const & conf,
		     simulation::Simulation const & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);
    
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);

    /**
     * store results in configuration
     */
    virtual int update_configuration(topology::Topology const & topo,
				     configuration::Configuration & conf,
				     simulation::Simulation const & sim);
    

    /**
     * calculate the interaction for a given atom pair.
     * SLOW! as it has to create the periodicity...
     */
    virtual int calculate_interaction(topology::Topology & topo,
			      configuration::Configuration & conf,
			      simulation::Simulation & sim,
			      unsigned int atom_i, unsigned int atom_j,
			      math::Vec & force, 
			      double &e_lj, double &e_crf);

    /**
     * calculate the hessian for a given atom.
     */
    virtual int calculate_hessian(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  unsigned int atom_i, unsigned int atom_j,
				  math::Matrix & hessian);

    Storage & shortrange_storage()
    {
      return m_storage;
    }
    Storage & longrange_storage()
    {
      return m_longrange_storage;
    }
    
  protected:
    Storage m_longrange_storage;

    Nonbonded_Outerloop m_outerloop;
  };
  
} // interaction

#endif
