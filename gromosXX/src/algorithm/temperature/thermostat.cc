/**
 * @file thermostat.cc
 * methods of the thermostat base class
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>
#include <configuration/state_properties.h>

#include "thermostat.h"

#include <util/debug.h>

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

      for(unsigned int i=last_atom + 1; i <= r_it->last_atom; ++i){
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
      unsigned int ir_bath, com_bath;

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
	  DEBUG(10, "com_v=" << math::v2s(new_com_v) << " ir_v=" << math::v2s(ir_v));

	} // loop over atoms
      } // loop over molecules of bath
    } // seperately coupled

    last_atom = r_it->last_atom;
    last_mol = r_it->last_molecule;
  } // loop over baths
}
