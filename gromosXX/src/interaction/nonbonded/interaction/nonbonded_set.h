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

namespace topology
{
  class Topology;
}
namespace configuration
{
  class Configuration;
}
namespace simulation
{
  class Simulation;
}

namespace interaction
{
  class Pairlist_Algorithm;
  class Nonbonded_Parameter;
  
  /**
   * @class Nonbonded_Set
   * calculates the nonbonded interactions.
   */
  class Nonbonded_Set
  {
  public:    
    /**
     * Constructor.
     */
    Nonbonded_Set(Pairlist_Algorithm & pairlist_alg, Nonbonded_Parameter & param);

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
		     bool quiet = false);
    
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim,
				       int tid = 0, int num_threads = 1);

    /**
     * calculate the interaction for a given atom pair.
     * SLOW! as it has to create the periodicity...
     */
    int calculate_interaction(topology::Topology & topo,
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
      return m_shortrange_storage;
    }
    Storage & longrange_storage()
    {
      return m_longrange_storage;
    }
    
    Pairlist & pairlist() { return m_pairlist; }
    Pairlist const & pairlist()const { return m_pairlist; }

  protected:
    Pairlist_Algorithm & m_pairlist_alg;
    
    Storage m_shortrange_storage;
    Storage m_longrange_storage;

    Pairlist m_pairlist;

    Nonbonded_Outerloop m_outerloop;

  };
  
} // interaction

#endif
