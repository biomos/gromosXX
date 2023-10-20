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
 * @file perturbed_cg_bond_interaction.h
 * perturbed cg bond interaction.
 */

#ifndef INCLUDED_PERTURBED_CG_BOND_INTERACTION
#define INCLUDED_PERTURBED_CG_BOND_INTERACTION

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
   * @class Perturbed_CG_Bond_Interaction
   * calculates the perturbed bond interactions (cg).
   */
  class Perturbed_DP_Bond_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_DP_Bond_Interaction(
       DP_Bond_Interaction &bond_interaction)
      : Interaction("PerturbedDPBond"),
	m_interaction(bond_interaction)
    {}
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_DP_Bond_Interaction() {}

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
      // os << "Perturbed cg bond interaction\n";
      return 0;
    };

    /**
     * calculate the interactions (force and energy, lambda derivative)
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
    
  protected:
    DP_Bond_Interaction & m_interaction;

  };
  
} // interaction

#endif
