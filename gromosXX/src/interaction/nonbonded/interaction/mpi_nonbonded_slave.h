/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

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
