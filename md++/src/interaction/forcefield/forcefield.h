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
  class Forcefield : public std::vector<Interaction *>,
		     public algorithm::Algorithm
  {
  public:
    /**
     * Constructor
     */
    Forcefield() 
      : std::vector<Interaction *>(),
	algorithm::Algorithm("Forcefield") {}
    /**
     * Destructor
     */
    ~Forcefield();
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

    /**
     * GPU-native Interactions (NonBonded, bonded terms) atomicAdd into
     * the mirror's shared force/virial buffers and leave them resident
     * (calculate_interactions() zeros them once via sim.cuda().
     * zero_mirror_force(), not per-Interaction) -- exempting FORCE/
     * VIRIAL here keeps Algorithm_Sequence::run()'s default post-
     * apply() invalidation from immediately erasing that freshness,
     * letting Leap_Frog_Velocity<gpuBackend> read the accumulated
     * force with zero extra round trip. Everything else (POS/VEL/BOX)
     * keeps the default MIRROR_ALL behaviour -- Forcefield never
     * writes those itself.
     */
    virtual unsigned gpu_mirror_touches() const override {
      return gpu::MIRROR_ALL & ~(gpu::MIRROR_FORCE | gpu::MIRROR_VIRIAL);
    }

  protected:

  };
  
} // interaction

#endif
