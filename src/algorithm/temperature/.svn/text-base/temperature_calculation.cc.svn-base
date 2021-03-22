/**
 * @file temperature_calculation.cc
 * calculates the temperature.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../configuration/state_properties.h"

#include "../../io/print_block.h"

#include "temperature_calculation.h"

#undef MODULE
#undef SUBMODULE

#define MODULE algorithm
#define SUBMODULE temperature

#include "../../util/debug.h"

int algorithm::Temperature_Calculation
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  DEBUG(7, "Temperature calculation");
  
  m_timer.start();
  
  // zero previous (temperature scaling) energies
  conf.old().energies.zero(false, true);
  // zero the energies in the multibath
  DEBUG(8, "\tbaths: " << unsigned(sim.multibath().size()));
  
  for(unsigned int i=0; i < unsigned(sim.multibath().size()); ++i)
    sim.multibath().bath(i).ekin = 0.0;

  topology::Temperaturegroup_Iterator tg_it = topo.temperature_group_begin(),
    tg_to = topo.temperature_group_end();
  
  math::Vec com_v, new_com_v;
  double com_ekin, ekin, new_com_ekin, new_ekin;
  
  unsigned int ir_bath, com_bath;
  
  configuration::State_Properties state_props(conf);

  for( ; tg_it != tg_to; ++tg_it){
    state_props.
      molecular_translational_ekin(sim, tg_it.begin(), tg_it.end(),
				   topo.mass(),
				   com_v, com_ekin, ekin,
				   new_com_v, new_com_ekin, new_ekin);

    DEBUG(9, "temperature group: " << *tg_it.begin() << " - " << *tg_it.end());
    DEBUG(10, "average com_v:" << math::v2s(com_v) << " com_ekin:" 
	  << com_ekin << " ekin:" << ekin);
    
    DEBUG(10, "new com_v:" << math::v2s(new_com_v) << " com_ekin:" 
	  << new_com_ekin << " ekin:" << new_ekin);
    
    sim.multibath().in_bath(*tg_it.begin(), com_bath, ir_bath);

    DEBUG(15, "adding to bath: com: "
	  << com_bath << " ir: " << ir_bath);
    DEBUG(20, "number of baths: energy " 
	  << unsigned(conf.old().energies.kinetic_energy.size())
	  << " com ekin "
	  << unsigned(conf.old().energies.com_kinetic_energy.size())
	  << " ir ekin "
	  << unsigned(conf.old().energies.ir_kinetic_energy.size())
	  );
    
    // store the new ones in multibath (velocity scaling for next step)
    sim.multibath().bath(com_bath).ekin += new_com_ekin;
    sim.multibath().bath(ir_bath).ekin += new_ekin - new_com_ekin;

    // and the averages in the energies
    conf.old().energies.com_kinetic_energy[com_bath] += com_ekin;
    conf.old().energies.ir_kinetic_energy[ir_bath] += ekin - com_ekin;

  }

  // loop over the bath kinetic energies
  for(size_t i=0; i<conf.old().energies.kinetic_energy.size(); ++i)
    conf.old().energies.kinetic_energy[i] =
      conf.old().energies.com_kinetic_energy[i] +
      conf.old().energies.ir_kinetic_energy[i];

  // and the perturbed energy derivatives (if there are any)
  if (sim.param().perturbation.perturbation){
    
    math::VArray &vel = conf.current().vel;
    math::VArray &old_vel = conf.old().vel;
  
    // loop over the baths
    std::vector<simulation::bath_index_struct>::iterator
      it = sim.multibath().bath_index().begin(),
      to = sim.multibath().bath_index().end();
  
    unsigned int last = 0;
    std::vector<double> &e_kin = 
      conf.old().perturbed_energy_derivatives.kinetic_energy;

    assert(e_kin.size() == conf.old().energies.kinetic_energy.size());

    // zero the kinetic energy derivatives
    for(unsigned int i=0; i< e_kin.size(); i++){
      e_kin[i] = 0.0;
    }

    for(; it != to; ++it){
      DEBUG(10, "starting atom range...");
      
      // or just put everything into the first bath...?
      unsigned int bath = it->com_bath;
    
      // don't zero them here, if there is a second set of DOF coming that is coupled to the 
      // same bath, you loose your contribution so far.
      //e_kin[bath] = 0.0;
      DEBUG(10, "\tperturbed kinetic energy last atom: "
	    << it->last_atom << " to bath " << bath << "\n");
      
      for(unsigned int i=last; i<=it->last_atom; ++i){
	
	if (i >= topo.num_solute_atoms()) break;
	
	if(topo.is_perturbed(i)){
	  DEBUG(11, "\tperturbed kinetic energy for " << i 
		<< " in bath " << bath);
	  DEBUG(11, "\tA_mass: " << topo.perturbed_solute().atoms()[i].A_mass() 
		<< " B_mass: " << topo.perturbed_solute().atoms()[i].B_mass());
	  DEBUG(11, "\tabs2(vel): " << math::abs2(vel(i)));
	  DEBUG(11, "\tabs2(old_vel): " << math::abs2(old_vel(i)));
	  
	  DEBUG(10, "\tbefore dE_kin/dl: " << e_kin[bath]);

	  const double l_deriv = 
	    topo.individual_lambda_derivative(simulation::mass_lambda)
	    [topo.atom_energy_group()[i]][topo.atom_energy_group()[i]];
	
          DEBUG(10, "\tmass:    " << topo.mass()[i]);
	  
	  DEBUG(10, "\tl_deriv: " << l_deriv);
	  
	  // for some reason we take the new velocities here
          // let's take the average of the old and the new velocities (squared) to be consistent with the 
          // calculation of the kinetic energy itself
	  e_kin[bath] -= l_deriv * 
	    (topo.perturbed_solute().atoms()[i].B_mass() -
	     topo.perturbed_solute().atoms()[i].A_mass()) 
	    * (math::abs2(old_vel(i)) + math::abs2(vel(i))) * 0.5;
	  
	  DEBUG(10, "\tdE_kin/dl: " << e_kin[bath]);
	
	}
      
      } // atoms in bath
    
      last = it->last_atom + 1;
    
    } // bath indices!!!!

    for(size_t i=0; i < e_kin.size(); ++i){
      e_kin[i] *= 0.5;
    }
    
  } // if perturbation...

  m_timer.stop();
  
  return 0;
  
}

int algorithm::Temperature_Calculation
::init(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation & sim,
       std::ostream & os,
       bool quiet)
{
  apply(topo, conf, sim);
  
  if (!quiet){
    // os << "Temperature calculation\n";
    io::print_MULTIBATH_COUPLING(os, sim.multibath());
    io::print_DEGREESOFFREEDOM(os, sim.multibath());
    io::print_MULTIBATH(os, sim.multibath(),
			conf.old().energies,
			"INITIAL TEMPERATURES");
  }
  return 0;
}

