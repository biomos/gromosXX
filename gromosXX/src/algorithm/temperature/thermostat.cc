/**
 * @file thermostat.cc
 * methods of the thermostat base class
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../configuration/state_properties.h"

#include "thermostat.h"

#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE

#define MODULE algorithm
#define SUBMODULE temperature

void algorithm::Thermostat
::scale(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  math::VArray &vel = conf.current().vel;
  
  // loop over the ranges
  std::vector<simulation::bath_index_struct>::const_iterator
    r_it = sim.multibath().bath_index().begin(),
    r_to = sim.multibath().bath_index().end();

  int last_atom = -1, last_tg = -1;
  configuration::State_Properties state_props(conf);

  for(; r_it != r_to; ++r_it){

    DEBUG(8, "atoms " << last_atom + 1 << " - " << r_it->last_atom
	  << " temperature groups " << last_tg + 1 << " - " << r_it->last_temperature_group);
    
    // decide whether molecular translational kinetic energy and
    // internal, rotational kinetic energies are jointly coupled
    if (r_it->com_bath == r_it->ir_bath){
      if (sim.multibath()[r_it->com_bath].scale == 1.0) continue;
      
      DEBUG(8, "jointly coupled, scaling with "
	    << sim.multibath()[r_it->com_bath].scale);

      for(unsigned int i=last_atom + 1; i <= r_it->last_atom; ++i){
	vel(i) *= sim.multibath()[r_it->com_bath].scale;
      }

    }
    else{
  
      if (sim.multibath()[r_it->com_bath].scale == 1.0 && 
	  sim.multibath()[r_it->ir_bath].scale == 1.0) continue;

      topology::Temperaturegroup_Iterator 
	tg_it = topo.temperature_group_begin(),
	tg_to = topo.temperature_group_begin();
      
      tg_it += last_tg + 1;
      tg_to += r_it->last_temperature_group + 1;

      DEBUG(8, "separately coupled");
      DEBUG(8, "scaling from temperature group " << last_tg << " to "
	    << r_it->last_temperature_group);

      math::Vec com_v, new_com_v, ir_v;
      double com_ekin = 0.0, ekin = 0.0, new_com_ekin = 0.0, new_ekin = 0.0;
      unsigned int ir_bath = 0, com_bath = 0;

      // which bathes?
      sim.multibath().in_bath(*(tg_it.begin()), com_bath, ir_bath);  

      for( ; tg_it != tg_to; ++tg_it){
	// new molecular translational velocities
	state_props.
	  molecular_translational_ekin(sim, tg_it.begin(), tg_it.end(),
				       topo.mass(),
				       com_v, com_ekin, ekin,
				       new_com_v, new_com_ekin, new_ekin);

	topology::Atom_Iterator start = tg_it.begin(),
	  end = tg_it.end();
    
	// loop over the atoms in the temperature group
	for( ; start != end; ++start){
    
	  assert(unsigned(vel.size()) > *start);
	  ir_v = vel(*start) - new_com_v;
	  vel(*start) =
	    sim.multibath()[com_bath].scale * new_com_v +
	    sim.multibath()[ir_bath].scale * ir_v;
	  
	  DEBUG(10, "com scale=" << sim.multibath()[com_bath].scale
		<< " ir scale=" << sim.multibath()[ir_bath].scale);
	  DEBUG(10, "com_v=" << math::v2s(new_com_v) << " ir_v=" << math::v2s(ir_v));

	} // loop over atoms
      } // loop over temperature group of bath
    } // seperately coupled

    last_atom = r_it->last_atom;
    last_tg = r_it->last_temperature_group;
  } // loop over baths
}

