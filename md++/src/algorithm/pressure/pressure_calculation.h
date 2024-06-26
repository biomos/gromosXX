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
 * @file pressure_calculation.h
 * pressure calculation
 */

#ifndef INCLUDED_PRESSURE_CALCULATION_H
#define INCLUDED_PRESSURE_CALCULATION_H

namespace algorithm
{
  
  /**
   * @class Pressure_Calculation
   * pressure calculation.
   */
  class Pressure_Calculation : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Pressure_Calculation() : Algorithm("PressureCalculation") {}

    /**
     * Destructor.
     */
    virtual ~Pressure_Calculation() {}
    
    /**
     * apply the pressure calculation
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim);
        
    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false) 
    {
      // os << "Pressure calculation\n";
      return 0;
    };

  private:

  };
  
} // algorithm

#endif
