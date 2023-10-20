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
 * @file leap_frog.h
 * the leap frog algorithm.
 */

#ifndef INCLUDED_LEAP_FROG_H
#define INCLUDED_LEAP_FROG_H

namespace algorithm
{
  /**
   * @class Leap_Frog_Position
   * implements the leap frog algorithm for the positions.
   */
  class Leap_Frog_Position : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Leap_Frog_Position() : Algorithm("Leap_Frog_Position") {}

    /**
     * Destructor.
     */
    virtual ~Leap_Frog_Position(){}
    
    /**
     * Leap frog step.
     */
    virtual int apply(topology::Topology &topo, 
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false) 
    {
      if (!quiet)
	os << "\tLeap frog position\nEND\n";
      return 0;
    };

  };

  /**
   * @class Leap_Frog_Velocity
   * implements the leap frog algorithm for the velocities.
   */
  class Leap_Frog_Velocity : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Leap_Frog_Velocity() : Algorithm("Leap_Frog_Velocity") {};
    /**
     * Destructor.
     */
    virtual ~Leap_Frog_Velocity(){};
    
    /**
     * Leap frog step.
     */
    virtual int apply(topology::Topology &topo, 
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false) 
    {
      if (!quiet)
	os << "INTEGRATION\n\tLeap frog velocity\n";
      return 0;
    };

  };

  
} // algorithm

#endif

