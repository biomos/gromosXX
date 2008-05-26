/**
 * @file berendsen_barostat.cc
 * methods of the berendsen barostat.
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <configuration/state_properties.h>

#include <util/error.h>

#include "berendsen_barostat.h"

#undef MODULE
#undef SUBMODULE

#define MODULE algorithm
#define SUBMODULE pressure

#include <util/debug.h>

int algorithm::Berendsen_Barostat
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  m_timer.start();
  
  DEBUG(8, "Berendsen Barostat == apply");

  // position are current!
  math::VArray & pos = conf.current().pos;
  math::Matrix & pressure = conf.old().pressure_tensor;
  math::Box & box = conf.current().box;

  DEBUG(9, "scaling = " << sim.param().pcouple.scale);
  
  switch(sim.param().pcouple.scale){
    case math::pcouple_isotropic:
      {
	DEBUG(9, "total pressure (isotropic) ...");
	double total_pressure =  (pressure(0,0)
				  + pressure(1,1)
				  + pressure(2,2)) / 3.0;

	DEBUG(8, "pressure: " << total_pressure);
	
	double mu = pow(1.0 - sim.param().pcouple.compressibility
			* sim.time_step_size() / sim.param().pcouple.tau
			* (sim.param().pcouple.pres0(0,0) - total_pressure),
			1.0/3.0);

	DEBUG(8, "mu: " << mu);

	// scale the box
	box *= mu;

	// scale the positions
	for(unsigned int i=0; i<pos.size(); ++i)
	  pos(i) = mu * pos(i);

        // scale the reference positions
        if (sim.param().posrest.scale_reference_positions){
          std::vector<topology::position_restraint_struct>::iterator
              it = topo.position_restraints().begin(),
              to = topo.position_restraints().end();
          for(; it != to; ++it) 
            it->pos *= mu;
        }

	break;
      }
    case math::pcouple_anisotropic:
      {
	math::Vec mu;

	DEBUG(8, "anisotropic pressure scaling");

	for(int i=0; i<3; ++i){
	  mu(i) = pow(1.0 - sim.param().pcouple.compressibility
		      * sim.time_step_size() / sim.param().pcouple.tau
		      * (sim.param().pcouple.pres0(i,i) - 
			 pressure(i,i)),
		      1.0/3.0);
	}

	DEBUG(10, "mu = " << math::v2s(mu));
	
	// scale the box
	for(int i=0; i<3; ++i)
	  for(int j=0; j<3; ++j)
	    box(i)(j) *= mu(j);

	DEBUG(10, "and the positions...");

	// scale the positions
	for(unsigned int i=0; i<pos.size(); ++i)
	  for(int j=0; j<3; ++j)
	    pos(i)(j) *= mu(j);
        
	// scale the reference positions
        if (sim.param().posrest.scale_reference_positions){
          std::vector<topology::position_restraint_struct>::iterator
              it = topo.position_restraints().begin(),
              to = topo.position_restraints().end();
          for(; it != to; ++it) 
            for(int j=0; j<3; ++j)
              it->pos(j) *= mu(j);
        }
        
	break;
      }
    case math::pcouple_full_anisotropic:
      {
	
	math::Matrix mu;

	for(int i=0; i<3; ++i){
	  for(int j=0; j<3; ++j){
	  
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
	for(unsigned int i=0; i<pos.size(); ++i)
	  pos(i) = math::product(mu, pos(i));
        
        // scale the reference positions
        if (sim.param().posrest.scale_reference_positions){
          std::vector<topology::position_restraint_struct>::iterator
              it = topo.position_restraints().begin(),
              to = topo.position_restraints().end();
          for(; it != to; ++it) 
            it->pos = math::product(mu, it->pos);
        }

      }
    default:
      return 0;
  }

  // check periodicity
  switch(conf.boundary_type){
    case math::vacuum:
      break;
    case math::rectangular:
      {
        
	if (abs(conf.current().box(0)) <= 2*sim.param().pairlist.cutoff_long ||
	    abs(conf.current().box(1)) <= 2*sim.param().pairlist.cutoff_long ||
	    abs(conf.current().box(2)) <= 2*sim.param().pairlist.cutoff_long){
	  io::messages.add("box is too small: not twice the cutoff!",
			   "configuration",
			   io::message::error);
	}
	
	break;
      }
    case math::triclinic:
      {
	// NO CUTOFF CHECK -- IMPLEMENT!!!
	break;
      }
    case math::truncoct:
      {
	if (0.5 * sqrt(3.0) * abs(conf.current().box(0)) <= 2 * sim.param().pairlist.cutoff_long){
	  
	  io::messages.add("box is too small: not 4 / sqrt(3) * cutoff!",
			   "configuration",
			   io::message::critical);
	  return E_BOUNDARY_ERROR;
	}
	break;
      }
    default:
      std::cout << "wrong periodic boundary conditions!";
      io::messages.add("wrong PBC!", "In_Configuration", io::message::error);
  }

  m_timer.stop();

  return 0;
  
}
