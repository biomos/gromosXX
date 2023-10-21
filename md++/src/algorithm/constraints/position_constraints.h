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
 * @file position_constraints.h
 * position constraints
 */

#ifndef INCLUDED_POSITION_CONSTRAINTS_H
#define INCLUDED_POSITION_CONSTRAINTS_H

namespace algorithm
{
  /**
   * @class Position_Constraints
   * implements position constraints
   */
  class Position_Constraints : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Position_Constraints(std::string const name = "Position_Constraints");

    /**
     * Destructor.
     */
    virtual ~Position_Constraints();
    
    /**
     * apply position constraints
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim);

    /**
     * initialise
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);

  protected:

  };
  
} //algorithm

#endif


