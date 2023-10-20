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
 * @file leap_frog.cc
 * contains the implementation
 * for the classes Leap_Frog_Position and Leap_Frog_Velocity.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "leap_frog.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE integration

/**
 * Leap frog step.
 */
int algorithm::Leap_Frog_Position
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation &sim)
{
  m_timer.start(sim);
  
  const int num_atoms = topo.num_atoms();
  
  // r = r + v*dt
#ifdef OMP
#pragma omp parallel for
#endif
  for(int i=0; i < num_atoms; ++i)
    conf.current().pos(i) =
      conf.old().pos(i) + conf.current().vel(i) * sim.time_step_size();
  
  if (sim.param().polarise.cos) {
#ifdef OMP
#pragma omp parallel for
#endif
    for(int i=0; i < num_atoms; ++i) {
      // Verlet type prediction of cos displacement 
      // conf.current() contains the oldold posV
      conf.current().posV(i) = 2*conf.old().posV(i) - conf.current().posV(i);

    }
  }
  
  m_timer.stop();
  
  return 0;
}

/**
 * Leap frog step.
 */
int algorithm::Leap_Frog_Velocity
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
#ifdef OMP
#pragma omp parallel for
#endif
  for(int i=0; i < num_atoms; ++i){
    conf.current().vel(i) =
      conf.old().vel(i) + conf.old().force(i) * sim.time_step_size() / topo.mass()(i);

    DEBUG(10, "atom " << i
	  << "\n\tf=" << math::v2s(conf.old().force(i))
	  << "\n\tmass=" << topo.mass()(i)
	  << "\n\tvel=" << math::v2s(conf.old().vel(i)));
  }
  
  m_timer.stop();
  
  return 0;
  
}
