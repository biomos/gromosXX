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
