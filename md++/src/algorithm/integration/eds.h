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

/* 
 * File:   eds.h
 * Author: haniels
 *
 * Created on August 2, 2011, 1:41 PM
 */

#ifndef EDS_H
#define	EDS_H

#include "H_acceleration.h"
#include <numeric>
namespace algorithm
{
  /**
   * @class EDS
   * implements EDS.
   */
  class EDS : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    EDS() : Algorithm("EDS"), conf2(NULL) {}
    
    void set_conf2(configuration::Configuration & conf) {
      conf2 = &conf;
    }
    

    /**
     * Destructor.
     */
    virtual ~EDS(){}
    
    /**
     * 
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
		     bool quiet = false);
    
    /**
     * Check for round trips
     **/
    bool check_round_trip(simulation::Simulation &sim);
    
    /**
     * get average over offsets
     **/
    double getAverage(simulation::Simulation &sim);

    AccelerationContainer accel_cont; // stores acceleration parameters
    EnergyAcceleration * accel_ptr; // stores acceleration parameters (for single-site) -> will be depricated
    
   private:
     configuration::Configuration * conf2;
  
  };
   
} // algorithm

#endif



