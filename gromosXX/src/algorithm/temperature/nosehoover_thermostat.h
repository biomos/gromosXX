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
 * @file nosehoover_thermostat.h
 * Nose-Hoover Thermostat.
 */

#ifndef INCLUDED_TEMPERATURE_NOSEHOOVER_H
#define INCLUDED_TEMPERATURE_NOSEHOOVER_H

#include "thermostat.h"

namespace algorithm
{
  
  /**
   * @class NoseHoover_Thermostat
   * the Nose-Hoover Thermostat.
   */
  class NoseHoover_Thermostat : public Thermostat
  {
  public:
    /**
     * Constructor.
     */
    NoseHoover_Thermostat()
      : Thermostat("BerendsenThermostat") {}

    /**
     * Destructor.
     */
    virtual ~NoseHoover_Thermostat() {}
    
    /**
     * apply the temperature scaling
     * for baths with tau=-1 nothing is done.
     * the kinetic energy can not be calculated here, because
     * later on SHAKE might be applied.
     * @param topo the Topology
     * @param conf the Configuration
     * @param sim the Simulation
     */
    virtual int apply
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim
     );

    /**
     * initialise
     */
    virtual int init
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     std::ostream & os = std::cout,
     bool quiet = false
     );

    /**
     * calculate the scaling factors.
     */
    void calc_scaling
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim
     );

    /**
     * calculate the scaling factors for Nose-Hoover Chains
     */
    void calc_chain_scaling
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim
     );

  private:

  };
  
} // algorithm

#endif
