/**
 * @file leap_frog.cc
 * contains the implementation
 * for the classes Leap_Frog_Position and Leap_Frog_Velocity.
 */

#include <util/stdheader.h>

#include <topology/core/core.h>
#include <topology/topology.h>
#include <simulation/multibath.h>
#include <simulation/parameter.h>
#include <simulation/simulation.h>
#include <configuration/energy.h>
#include <configuration/energy_average.h>
#include <configuration/configuration.h>

#include <algorithm/algorithm.h>
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
  // r = r + v*dt
  conf.current().pos = conf.old().pos + conf.current().vel
    * sim.time_step_size();
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
  conf.exchange_state();
  // copy the box
  conf.current().box = conf.old().box;
  
  // v = v + f * dt / m
  conf.current().vel = conf.old().vel + conf.old().force 
    * sim.time_step_size() / topo.mass();
  return 0;
  
}
