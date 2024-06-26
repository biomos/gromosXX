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
 * @file slow_growth.h
 * slow growth algorithm
 */

#ifndef INCLUDED_SLOW_GROWTH_H
#define INCLUDED_SLOW_GROWTH_H

namespace algorithm
{
  /**
   * @class Slow_Growth
   * implements slow growth.
   * mainly updates the topology for a changed
   * lambda value.
   */
  class Slow_Growth : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Slow_Growth() : Algorithm("SlowGrowth") {}

    /**
     * Destructor.
     */
    virtual ~Slow_Growth() {}
    
    /**
     * do the lambda change.
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
      os << "Slow growth\n";
      return 0;
    };

  };
  
} // algorithm

#endif
