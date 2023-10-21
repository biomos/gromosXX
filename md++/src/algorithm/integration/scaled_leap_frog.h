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
 * @file scaled_leap_frog.h
 * the leap frog algorithm.
 */

#ifndef INCLUDED_SCALED_LEAP_FROG_H
#define INCLUDED_SCALED_LEAP_FROG_H

namespace algorithm
{
  /**
   * @class Scaled_Leap_Frog_Velocity
   * implements the leap frog algorithm for the velocities.
   */
  class Scaled_Leap_Frog_Velocity : public Leap_Frog_Velocity
  {
  public:
/**
     * Constructor.
     */
    Scaled_Leap_Frog_Velocity() :
    Leap_Frog_Velocity() { name = "Scaled_Leap_Frog_Velocity"; }
    /**
     * Destructor.
     */
    virtual ~Scaled_Leap_Frog_Velocity(){};
    
    /**
     * Scaled leap frog step.
     */
    virtual int apply(topology::Topology &topo, 
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);


  };

  
} // algorithm

#endif

