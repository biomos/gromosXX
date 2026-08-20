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
 * @file energy_calculation.h
 * calculates the (total) energies and updates the averages
 */

#ifndef INCLUDED_ENERGY_CALCULATION_H
#define INCLUDED_ENERGY_CALCULATION_H

namespace algorithm
{
  /**
   * @class Energy_Calculation
   * calculates total energies, updates the averages
   */
  class Energy_Calculation : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Energy_Calculation() : Algorithm("EnergyCalculation") {}

    /**
     * Destructor.
     */
    virtual ~Energy_Calculation(){}
    
    /**
     * calculate the totals
     * update the averages
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

    // conf.old().energies.calculate_totals() (in apply()) reads
    // kinetic_energy, which Temperature_Calculation<gpuBackend> may
    // have only queued asynchronously on the GPU rather than finished
    // computing on the host yet -- see Algorithm::needs_finalized_gpu_
    // state()'s doc comment.
    virtual bool needs_finalized_gpu_state() const override { return true; }

  };
  
} // algorithm

#endif

