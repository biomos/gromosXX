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
 * @file perturbed_dihedral_new_interaction.h
 * perturbed  dihedral interaction.
 */

#ifndef INCLUDED_PERTURBED_DIHEDRAL_NEW_INTERACTION
#define INCLUDED_PERTURBED_DIHEDRAL_NEW_INTERACTION

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
   * @class Perturbed_Dihedral_new_Interaction
   * calculates the perturbed dihedral interactions.
   */
  class Perturbed_Dihedral_new_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Dihedral_new_Interaction(Dihedral_new_Interaction & dihedral_new_interaction)
      : Interaction("PerturbedDihedral"),
	m_interaction(dihedral_new_interaction)
    {}
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Dihedral_new_Interaction() {}

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
      // os << "Perturbed dihedral interaction\n";
      return 0;
    };
    /**
     * calculate the interactions (force and energy, lambda derivative)
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
    
  protected:
    Dihedral_new_Interaction & m_interaction;
  };
  
} // interaction

#endif
