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
 * @file forcefield.h
 * the Forcefield class.
 */

#ifndef INCLUDED_FORCEFIELD_H
#define INCLUDED_FORCEFIELD_H

namespace configuration{
	class Configuration;
}
namespace topology{
	class Topology;
}
namespace simulation{
	class Simulation;
}

/**
 * @namespace interaction
 * namespace that contains the classes to
 * handle the interactions between the particles.
 * (energies, forces).
 */
namespace interaction
{
  /**
   * @class Forcefield
   * contains the specific interactions.
   * (clear does not call them (i guess) -- sorry, don't know what this means anymore)
   * Strategy Pattern.
   */
  template<typename Backend = util::cpuBackend>
  class ForcefieldT : public std::vector<Interaction *>,
		     public algorithm::AlgorithmT<Backend>
  {
  public:
    /**
     * Constructor
     */
    ForcefieldT() 
      : std::vector<Interaction *>(),
	algorithm::AlgorithmT<Backend>("Forcefield") {}
    /**
     * Destructor
     */
    ~ForcefieldT();
    /**
     * initialise
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);

    /**
     * calculate all interactions.
     */
    int calculate_interactions(topology::Topology & topo,
			       configuration::Configuration & conf,
			       simulation::Simulation & sim);


    /**
     * let the forcefield be used as an algorithm
     */
    int apply(topology::Topology & topo,
	      configuration::Configuration & conf,
	      simulation::Simulation & sim)
    {
      return calculate_interactions(topo, conf, sim);
    }
    
    virtual void print_timing(std::ostream & os);

    /**
     * const interaction accessor.
     */
    Interaction const * interaction(std::string name)const;

    /**
     * interaction accessor.
     */
    Interaction * interaction(std::string name);
    
  protected:

  };

  /**
   * @brief Allow use of Forcefield directly - defaults to ForcefieldT<util::cpuBackend>
   * 
   */
  using Forcefield = ForcefieldT<util::cpuBackend>;
  
} // interaction

#endif
