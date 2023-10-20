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
 * @file perturbed_soft_angle_interaction.h
 * perturbed soft angle interaction.
 */

#ifndef INCLUDED_PERTURBED_SOFT_ANGLE_INTERACTION
#define INCLUDED_PERTURBED_SOFT_ANGLE_INTERACTION

namespace configuration{
	class Configuration;
}
namespace topology{
	class Topology;
}
namespace simulation{
	class Simulation;
}

namespace interaction
{
  /**
   * @class Perturbed_Angle_Interaction
   * calculates the perturbed angle interactions.
   */
  class Perturbed_Soft_Angle_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Soft_Angle_Interaction()
      : Interaction("PerturbedSoftAngle") {}
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Soft_Angle_Interaction() {}

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false) ;
    /**
     * calculate the interactions (force and energy, lambda derivative)
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
    
  };
  
} // interaction

#endif
