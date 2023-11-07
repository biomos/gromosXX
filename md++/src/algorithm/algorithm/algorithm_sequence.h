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
 * @file algorithm_sequence.h
 * the Algorithm class.
 */

#ifndef INCLUDED_ALGORITHM_ALGORITHM_H
#define INCLUDED_ALGORITHM_ALGORITHM_H

namespace algorithm
{
  /**
   * @class Algorithm_Sequence
   * contains the specific algorithms.
   * An almost clean implementation of the Strategy pattern.
   */
  class Algorithm_Sequence : public std::vector<Algorithm *>
  {
  public:
    /**
     * Constructor
     */
    Algorithm_Sequence(bool clean = true);

    /**
     * Destructor
     */
    ~Algorithm_Sequence();

    /**
     * init
     */
    int init(topology::Topology &topo, 
	     configuration::Configuration &conf,
	     simulation::Simulation &sim,
	     std::ostream & os = std::cout,
	     bool quiet = false);

    /**
     * calculate all interactions.
     */
    int run(topology::Topology &topo, 
	    configuration::Configuration &conf,
	    simulation::Simulation &sim);

    /**
     * print timing information
     */
    int print_timing(std::ostream & os);

    /**
     * algorithm accessor
     */
    Algorithm * algorithm(std::string name);
    
    /**
     * print algorithm sequence.
     * this is a nice debugging function.
     */
    void printSequence();
    
  protected:
    bool clean;

  };
  
} // algorithm


#endif
