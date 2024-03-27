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
