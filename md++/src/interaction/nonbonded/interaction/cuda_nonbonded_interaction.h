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
 * @file cuda_nonbonded_interaction.h
 * the non bonded interactions:
 * Lennard-Jones and Coulomb interactions
 * using CUDA acceleration
 */

#pragma once

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

#include "interaction.h"
#include "nonbonded_parameter.h"

namespace interaction
{

  class Pairlist_Algorithm;
  class Nonbonded_Set_Interface;
  
  /**
   * @class CUDA_Nonbonded_Interaction
   * calculates the nonbonded interactions.
   */
  class CUDA_Nonbonded_Interaction : 
    public Nonbonded_Interaction
  {
  public:    

    /**
     * Constructor.
     */
    CUDA_Nonbonded_Interaction(Pairlist_Algorithm *pa);
    
    /**
     * Destructor.
     */
    virtual ~CUDA_Nonbonded_Interaction();

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

  };
  
} // interaction
