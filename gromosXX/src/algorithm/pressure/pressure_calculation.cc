/**
 * @file pressure_calculation.cc
 * calculates the pressure.
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>
#include <configuration/state_properties.h>

#include <math/volume.h>
#include <math/transformation.h>

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

  m_timer.start();
  // virial -0.5 factor applied here!
  conf.old().virial_tensor *= -0.5;
  // calculate the pressure tensor
  conf.old().pressure_tensor = topo.tot_cg_factor() *
          (conf.old().kinetic_energy_tensor - conf.old().virial_tensor) *
          (2.0 / math::volume(conf.old().box, conf.boundary_type));
  
  m_timer.stop();
  
  return 0;
  
}
