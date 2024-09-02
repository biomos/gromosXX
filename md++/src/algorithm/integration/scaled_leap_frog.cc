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
 * @file scaled_leap_frog.cc
 * contains the implementation
 * for the class  Scaled_Leap_Frog_Velocity.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "leap_frog.h"
#include "scaled_leap_frog.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE integration

/**
 * Leap frog step.
 */
int algorithm::Scaled_Leap_Frog_Velocity
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation &sim)
{
  m_timer.start(sim);

  conf.exchange_state();
  // copy the box
  conf.current().box = conf.old().box;

  const int num_atoms = topo.num_atoms();

  // v = v + f * dt / m
  int addc_index = 0;
  double scale_force = 0.0;
#ifdef OMP
#pragma omp parallel for
#endif
  for(int i=0; i < num_atoms; ++i){
  addc_index = sim.param().addecouple.check_index_adc(i);
    scale_force=1;
    if(addc_index!=-1)
      scale_force= sim.param().addecouple.adc_index()[addc_index].sv;
    
    conf.current().vel(i) =
      conf.old().vel(i) + scale_force * conf.old().force(i) * sim.time_step_size() 
      / topo.mass()(i);

    DEBUG(10, "atom " << i
	  << "\n\tf=" << math::v2s(conf.old().force(i))
	  << "\n\tmass=" << topo.mass()(i)
	  << "\n\tvel=" << math::v2s(conf.old().vel(i)));
  }
  
  m_timer.stop();
  
  return 0;
  
}
