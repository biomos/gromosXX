/**
 * @file temperature/berendsen_thermostat.tcc
 * methods of the berendsen thermostat
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
#include "berendsen_thermostat.h"

#undef MODULE
#undef SUBMODULE

#define MODULE algorithm
#define SUBMODULE temperature

#include <util/debug.h>

int algorithm::Berendsen_Thermostat
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  
  math::VArray &vel = conf.current().vel;
  // configuration::Energy &e = conf.current().energies;
  
  // loop over the baths
  std::vector<simulation::bath_struct>::iterator
    b_it = sim.multibath().begin(),
    b_to = sim.multibath().end();
  
  for(size_t num=0; b_it != b_to; ++b_it, ++num){
    // do temperature coupling for that bath?
    if (b_it->tau != -1){

      DEBUG(7, "pre-scale ekin: " << b_it->ekin);

      double free_temp = 2 * 
	b_it->ekin / (b_it->dof * math::k_Boltzmann);

      // divide by zero measure...
      if (free_temp < math::epsilon) free_temp = b_it->temperature;

      b_it->scale = sqrt(1.0 + sim.time_step_size() / b_it->tau *
			 (b_it->temperature / free_temp - 1));

      DEBUG(8, "free T " << free_temp << " dof " << b_it->dof);
      DEBUG(8, "ref. T " << b_it->temperature);
      DEBUG(8, "bath " << num << " scaling: " << b_it->scale);
      
    } // scaling ?
    else
      b_it->scale = 1;
  } // loop over the baths

  //--------------------------------
  // now we have the scaling factors
  //--------------------------------

  // loop over the ranges
  std::vector<simulation::bath_index_struct>::const_iterator
    r_it = sim.multibath().bath_index().begin(),
    r_to = sim.multibath().bath_index().end();

  int last_atom = -1, last_mol = -1;
  configuration::State_Properties state_props(conf);

  for(; r_it != r_to; ++r_it){

    DEBUG(8, "atoms " << last_atom + 1 << " - " << r_it->last_atom
	  << " mols " << last_mol + 1 << " - " << r_it->last_molecule);
    
    // decide whether molecular translational kinetic energy and
    // internal, rotational kinetic energies are jointly coupled
    if (r_it->com_bath == r_it->ir_bath){
      if (sim.multibath()[r_it->com_bath].tau == -1) continue;
      
      DEBUG(8, "jointly coupled, scaling with "
	    << sim.multibath()[r_it->com_bath].scale);

      for(size_t i=last_atom + 1; i <= r_it->last_atom; ++i){
	vel(i) *= sim.multibath()[r_it->com_bath].scale;
      }

    }
    else{
  
      topology::Molecule_Iterator 
	m_it = topo.molecule_begin(),
	m_to = topo.molecule_begin();
      
      m_it += last_mol + 1;
      m_to += r_it->last_molecule + 1;

      DEBUG(8, "separately coupled");
      DEBUG(8, "scaling from mol " << last_mol << " to "
	    << r_it->last_molecule);

      math::Vec com_v, new_com_v, ir_v;
      double com_ekin, ekin, new_com_ekin, new_ekin;
      size_t ir_bath, com_bath;

      // which bathes?
      sim.multibath().in_bath(*(m_it.begin()), com_bath, ir_bath);  

      for( ; m_it != m_to; ++m_it){
	// new molecular translational velocities
	state_props.
	  molecular_translational_ekin(m_it.begin(), m_it.end(),
				       topo.mass(),
				       com_v, com_ekin, ekin,
				       new_com_v, new_com_ekin, new_ekin);

	topology::Atom_Iterator start = m_it.begin(),
	  end = m_it.end();
    
	// loop over the atoms in the molecule
	for( ; start != end; ++start){
    
	  assert(unsigned(vel.size()) > *start);
	  ir_v = vel(*start) - new_com_v;
	  vel(*start) =
	    sim.multibath()[com_bath].scale * new_com_v +
	    sim.multibath()[ir_bath].scale * ir_v;
	  
	  DEBUG(10, "com scale=" << sim.multibath()[com_bath].scale
		<< " ir scale=" << sim.multibath()[ir_bath].scale);
	  DEBUG(10, "com_v=" << new_com_v << " ir_v=" << ir_v);

	} // loop over atoms
      } // loop over molecules of bath
    } // seperately coupled

    last_atom = r_it->last_atom;
    last_mol = r_it->last_molecule;
  } // loop over baths

  return 0;
  
}