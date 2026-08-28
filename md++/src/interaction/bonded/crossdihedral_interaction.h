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
 * @file crossdihedral_interaction.h
 * crossdihedral interaction.
 */

#ifndef INCLUDED_CROSSDIHEDRAL_INTERACTION_H
#define INCLUDED_CROSSDIHEDRAL_INTERACTION_H

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
   * @class Crossdihedral_Interaction
   * calculates the crossdihedral interactions.
   */
  class Crossdihedral_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Crossdihedral_Interaction() : Interaction("Crossdihedral") {}
    /**
     * Destructor.
     */
    virtual ~Crossdihedral_Interaction() {}

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
      // os << "Crossdihedral interaction\n";
      // sim.param().force.crossdihedral defaults to 1 with no FORCE-
      // block field to turn it off (simulation/parameter.h), so this
      // Interaction is pushed into every Forcefield regardless of
      // whether the topology defines any crossdihedral terms at all --
      // most topologies (anything without CMAP-style terms, e.g. this
      // ubiquitin benchmark) have zero. calculate_interactions() then
      // does nothing (an empty-range loop), but the base class's
      // needs_fresh_cpu_force()==true default still made Forcefield::
      // calculate_interactions() flush_gpu_dirty(MIRROR_FORCE) -- a
      // real device-wide sync + full force-array D2H copy -- before
      // every single call, for zero benefit. Cache emptiness once here
      // (topology is static for a run) instead of checking every step.
      m_has_terms = !topo.solute().crossdihedrals().empty();
      return 0;
    };
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);

    // See init()'s doc comment: only request the CPU-fresh force
    // publish when there are actually terms to read/write force for.
    virtual bool needs_fresh_cpu_force() const override { return m_has_terms; }

  protected:
    bool m_has_terms = true;

    /**
     * calculate nearest minimum
     */
    // double _calculate_nearest_minimum(double phi, int m, double pd);
    
  };
  
} // interaction

#endif
