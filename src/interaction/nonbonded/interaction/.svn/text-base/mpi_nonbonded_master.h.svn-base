/**
 * @file mpi_nonbonded_master.h
 * the non bonded interactions:
 * Lennard-Jones and Coulomb interactions.
 * MPI master:
 * distribute positions, calculate forces, collect forces
 */

#ifndef INCLUDED_MPI_NONBONDED_MASTER_H
#define INCLUDED_MPI_NONBONDED_MASTER_H

#include "nonbonded_interaction.h"

namespace interaction
{
  /**
   * @class MPI_Nonbonded_Master
   * calculates the nonbonded interactions using MPI
   * Master: distribute positions, calc forces, collect forces
   */
  class MPI_Nonbonded_Master : 
    public Nonbonded_Interaction
  {
  public:    
    /**
     * Constructor.
     */
    MPI_Nonbonded_Master(Pairlist_Algorithm *pa);
    /**
     * Destructor.
     */
    virtual ~MPI_Nonbonded_Master();
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

    /**
     * storage for stuff
     */
    Storage m_storage;

  };
  
} // interaction

#endif
