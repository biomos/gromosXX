/**
 * @file pressure/berendsen_barostat.cc
 * methods of the berendsen barostat.
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

#include <algorithm/algorithm.h>
#include "berendsen_barostat.h"

#undef MODULE
#undef SUBMODULE

#define MODULE algorithm
#define SUBMODULE temperature

#include <util/debug.h>

int algorithm::Berendsen_Barostat
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{

  math::VArray &pos = conf.old().pos;
  math::Matrix &pressure = conf.old().pressure_tensor;
  math::Box box = conf.old().box;

  switch(sim.param().pcouple.scale){
    case math::pcouple_isotropic:
      {
	double total_pressure =  (pressure(0,0)
				  + pressure(1,1)
				  + pressure(2,2)) / 3.0;
	
	double mu = pow(1.0 - sim.param().pcouple.compressibility
			* sim.time_step_size() / sim.param().pcouple.tau
			* (sim.param().pcouple.pres0(0,0) - total_pressure),
			1.0/3.0);

	// scale the box
	box = mu * box;

	// scale the positions
	for(int i=0; i<pos.size(); ++i)
	  pos(i) = mu * pos(i);

	break;
      }
    case math::pcouple_anisotropic:
      {
	math::Vec mu;

	for(int i=0; i<3; ++i){
	  mu(i) = pow(1.0 - sim.param().pcouple.compressibility
		      * sim.time_step_size() / sim.param().pcouple.tau
		      * (sim.param().pcouple.pres0(i,i) - 
			 pressure(i,i)),
		      1.0/3.0);
	}

	// scale the box
	for(int i=0; i<3; ++i)
	  box(i) = box(i) * mu;

	// scale the positions
	for(int i=0; i<pos.size(); ++i)
	  pos(i) = mu * pos(i);
	
	break;
      }
    case math::pcouple_full_anisotropic:
      {
	
	math::Matrix mu;

	for(int i=0; i<3; ++i){
	  for(int j=0; j<3; ++i){
	  
	    mu(i, j) = pow(1.0 - sim.param().pcouple.compressibility
			   * sim.time_step_size() / sim.param().pcouple.tau
			   * (sim.param().pcouple.pres0(i,j) -
			      pressure(i,j)),
			   1.0/3.0);
	  }
	}

	// scale the box
	box = math::product(mu, box);
	
	// scale the positions
	for(int i=0; i<pos.size(); ++i)
	  pos(i) = math::product(mu, pos(i));

      }
    default:
      return 0;
  }
  return 0;
  
}
