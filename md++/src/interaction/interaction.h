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
 * @file src/interaction/interaction.h
 * the interaction interface.
 */

#ifndef INCLUDED_INTERACTION_H
#define INCLUDED_INTERACTION_H

namespace configuration{
  class Configuration;
}
namespace topology{
  class Topology;
}
namespace simulation{
  class Simulation;
}
namespace util {
  class Algorithm_Timer;
}

namespace interaction
{
  /**
   * @class Interaction
   * @interface Interaction
   * declares the interaction interface.
   */
  class Interaction
  {
  public:
    /**
     * Constructor.
     */
    Interaction(std::string name) : name(name), m_timer(name) {};
    /**
     * Destructor.
     */
    virtual ~Interaction(){};
    /**
     * the name of the interaction.
     * can be used to identify a special class.
     */
    std::string name;
    /**
     * initialise
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false) = 0;
    // { return 0; }
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim) = 0;

    /**
     * @brief Does this Interaction read the *already-accumulated*
     * conf.current().force (e.g. to derive a correction term) before
     * computing/adding its own? Default false -- the overwhelming
     * majority of Interactions (bonded terms, NonBonded, restraint/
     * special forces) only ever += their own contribution, so they
     * don't care whether the array currently holds a CPU-fresh value or
     * is lagging behind an as-yet-unpublished GPU-resident write.
     * Molecular_Virial_Interaction is the one exception (it needs the
     * true per-atom total to correct atomic virial to molecular virial)
     * -- Forcefield::calculate_interactions() checks this flag and
     * publishes the GPU mirror's force before calling such an
     * Interaction, since GPU-native force writers (CUDA_Angle_
     * Interaction, CUDA_Nonbonded_Interaction, etc.) leave force
     * GPU-resident (gpu::MIRROR_FORCE marked dirty via mark_gpu_dirty(),
     * not synced back every call) until something downstream actually
     * needs the CPU-side value.
     */
    virtual bool needs_fresh_cpu_force() const { return false; }

    /**
     * timing information.
     */
    virtual void print_timing(std::ostream & os)
    {
      m_timer.print(os);
    }

  protected:
    /**
     * store time used in algorithm.
     */
    util::Algorithm_Timer m_timer;

  };  
  
} // interaction

#endif
