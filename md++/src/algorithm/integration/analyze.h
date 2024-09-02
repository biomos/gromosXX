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
 * @file analyze.h
 * re-evaluate a trajectory
 */

#ifndef INCLUDED_ANALYZE_H
#define INCLUDED_ANALYZE_H

namespace algorithm
{
  /**
   * @class Analyze_Step
   * implements analyzation of a trajectory
   */
  class Analyze_Step : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Analyze_Step(std::string trajectory);

    /**
     * Destructor.
     */
    virtual ~Analyze_Step(){}
    
    /**
     * analyze step.
     */
    virtual int apply(topology::Topology &topo, 
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);

    /**
     * init an algorithm
     * print out input parameter, what it does...
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);
    
  private:
    /**
     * trajectory to analyze
     */
    io::In_Configuration m_trajectory;

    /**
     * trajectory file
     */
   io::igzstream m_trajectory_file;
    
  };
  
} // algorithm

#endif
