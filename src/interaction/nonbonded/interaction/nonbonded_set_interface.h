/**
 * @file nonbonded_set_interface.h
 * the non bonded interactions for a set of atoms:
 * Lennard-Jones and Coulomb interactions
 * perturbed and non-perturbed.
 */

#ifndef INCLUDED_NONBONDED_SET_INTERFACE_H
#define INCLUDED_NONBONDED_SET_INTERFACE_H

#include "pairlist.h"
// #include "storage.h"
// #include "nonbonded_outerloop.h"

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
namespace util 
{
  class Algorithm_Timer;
}

namespace interaction
{
  class Pairlist_Algorithm;
  class Nonbonded_Parameter;
  class Storage;
  
  /**
   * @class Nonbonded_Set_Interface
   */
  class Nonbonded_Set_Interface
  {
  public:    
    /**
     * Constructor.
     */
    Nonbonded_Set_Interface(Pairlist_Algorithm & pairlist_alg, Nonbonded_Parameter & param,
			    int rank, int num_threads)
      : m_pairlist_alg(pairlist_alg), 
	m_rank(rank),
	m_num_threads(num_threads)
    {}

    /**
     * Destructor
     */
    virtual ~Nonbonded_Set_Interface() {}
    
    /**
     * initialize some things
     */
    virtual int init(topology::Topology const & topo,
		     configuration::Configuration const & conf,
		     simulation::Simulation const & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false) = 0;
    
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim) = 0;

    virtual int update_configuration(topology::Topology const & topo,
				     configuration::Configuration & conf,
				     simulation::Simulation const & sim) = 0;

    /**
     * calculate the interaction for a given atom pair.
     * SLOW! as it has to create the periodicity...
     */
    virtual int calculate_interaction(topology::Topology & topo,
				      configuration::Configuration & conf,
				      simulation::Simulation & sim,
				      unsigned int atom_i, unsigned int atom_j,
				      math::Vec & force, 
				      double &e_lj, double &e_crf) = 0;
    
    /**
     * calculate the hessian for a given atom.
     */
    virtual int calculate_hessian(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  unsigned int atom_i, unsigned int atom_j,
				  math::Matrix & hessian) = 0;

    Storage & storage(){ return m_storage; }

    PairlistContainer & pairlist() { return m_pairlist; }
    PairlistContainer const & pairlist()const { return m_pairlist; }

    void start_timer(std::string t) {
      if (m_rank == 0) m_pairlist_alg.timer().start(t);
    }
    void stop_timer(std::string t) {
      if (m_rank == 0) m_pairlist_alg.timer().stop(t);
    }

  protected:
    Pairlist_Algorithm & m_pairlist_alg;

    Storage m_storage;
    Storage m_storage_cuda;
    PairlistContainer m_pairlist;

    /**
     * OpenMP / MPI rank of thread running this set
     */
    int m_rank;
    /**
     * total number of OpenMP / MPI threads
     */
    int m_num_threads;
  };
  
} // interaction

#endif
