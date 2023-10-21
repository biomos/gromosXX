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
 * @file algorithm.h
 * base class for algorithms
 */

#ifndef INCLUDED_ALGORITHM_H
#define INCLUDED_ALGORITHM_H

namespace configuration
{
  class Configuration;
}
namespace topology
{
  class Topology;
}
namespace simulation
{
  class Simulation;
}
namespace util
{
  class Algorithm_Timer;
}

namespace algorithm
{
  /**
   * @class Algorithm
   * base class
   */
  class Algorithm
  {
  public:
    /**
     * Constructor.
     * @param name of the algorithm.
     */
    Algorithm(std::string name) : name(name), m_timer(name) {}

    /**
     * Destructor.
     */
    virtual ~Algorithm() {}
    
    /**
     * init an algorithm
     * print out input parameter, what it does...
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false) = 0;
    // { return 0; }
    
    /**
     * apply the algorithm
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim) {return 0;}

    /**
     * name of the algorithm
     */
    std::string name;

    /**
     * timing information.
     */
    virtual void print_timing(std::ostream & os)
    {
      m_timer.print(os);
    }
    /**
     * const accessor to timer
     */
    const util::Algorithm_Timer & timer() const {
      return m_timer;
    }
    /**
     * accessor to timer
     */
    util::Algorithm_Timer & timer() {
      return m_timer;
    }
    /**
     * accessor to timer
     */
    void timer(util::Algorithm_Timer &t) {
      m_timer = t;
    }
    
  protected:
    /**
     * store time used in algorithm.
     */
    util::Algorithm_Timer m_timer;
  };
}

#endif


