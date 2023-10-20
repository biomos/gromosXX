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
 * @file dihedral_new_interaction.h
 * dihedral interaction.
 */

#ifndef INCLUDED_DIHEDRAL_NEW_INTERACTION_H
#define INCLUDED_DIHEDRAL_NEW_INTERACTION_H

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
   * @class Dihedral_new_Interaction
   * calculates the dihedral interactions.
   */
  class Dihedral_new_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Dihedral_new_Interaction() : Interaction("Dihedral") {}
    /**
     * Destructor.
     */
    virtual ~Dihedral_new_Interaction() {}

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false) 
    {
      // if (!quiet)
      // os << "Dihedral interaction\n";
      return 0;
    };
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);

    
  protected:

    /**
     * calculate nearest minimum
     */
    // double _calculate_nearest_minimum(double phi, int m, double pd);
    
  };
  
} // interaction

#endif
