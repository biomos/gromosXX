/**
 * @file pressure_calculation.cc
 * calculates the pressure.
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
#include <configuration/state_properties.h>

#include <math/volume.h>

#include <algorithm/algorithm.h>
#include "pressure_calculation.h"

#undef MODULE
#undef SUBMODULE

#define MODULE algorithm
#define SUBMODULE temperature

#include <util/debug.h>

int algorithm::Pressure_Calculation
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  DEBUG(7, "Pressure calculation");

  // calculate the pressure tensor
  for(int i=0; i<3; ++i){
    for(int j=0; j<3; ++j){
      // virial -0.5 factor applied here!
      conf.old().virial_tensor(i,j) *= -0.5;
      conf.old().pressure_tensor(i,j) = 2 * 
	(conf.old().kinetic_energy_tensor(i,j)
	 - conf.old().virial_tensor(i,j)) /
	math::volume(conf.old().box, conf.boundary_type);
    }
  }
  
  // update the averages
  conf.current().energy_averages.update(conf.old().pressure_tensor,
					conf.old().energy_averages,
					sim.time_step_size());

  return 0;
  
}
