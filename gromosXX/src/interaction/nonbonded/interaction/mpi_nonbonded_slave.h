/**
 * @file mpi_nonbonded_slave.h
 * the non bonded interactions:
 * Lennard-Jones and Coulomb interactions.
 * MPI slave:
 * get positions, calculate forces and send them back
 */

#ifndef INCLUDED_MPI_NONBONDED_SLAVE_H
#define INCLUDED_MPI_NONBONDED_SLAVE_H

#include "nonbonded_interaction.h"

namespace interaction
{
  /**
   * @class MPI_Nonbonded_Slave
   * calculates the nonbonded interactions using MPI
   * Slave: receive positions, calc forces and send back to master
   */
  class MPI_Nonbonded_Slave : 
    public Nonbonded_Interaction
  {
  public:    
    /**
     * Constructor.
     */
    MPI_Nonbonded_Slave(Pairlist_Algorithm *pa);
    /**
     * Destructor.
     */
    virtual ~MPI_Nonbonded_Slave();
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
