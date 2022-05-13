/**
 * @file berendsen_thermostat.cc
 * methods of the berendsen thermostat
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../configuration/state_properties.h"

#include "berendsen_thermostat.h"

#undef MODULE
#undef SUBMODULE

#define MODULE algorithm
#define SUBMODULE temperature

#include "../../util/debug.h"

int algorithm::Berendsen_Thermostat::init
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 std::ostream & os,
 bool quiet
 )
{
  if (!quiet){
    os << "\tWeak-Coupling temperature coupling\n";
  }
  return 0;
}


int algorithm::Berendsen_Thermostat
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  m_timer.start();

  calc_scaling(topo, conf, sim);

  //--------------------------------
  // now we have the scaling factors
  //--------------------------------

  scale(topo, conf, sim);

  m_timer.stop();

  return 0;
  
}

void algorithm::Berendsen_Thermostat
::calc_scaling(topology::Topology & topo,
	       configuration::Configuration & conf,
	       simulation::Simulation & sim,
	       bool immediate)
{
  DEBUG(7, "Temperature Scaling!!!!");

  // std::cout.precision(9);
  // std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
  
  // loop over the baths
  std::vector<simulation::bath_struct>::iterator
    b_it = sim.multibath().begin(),
    b_to = sim.multibath().end();
  
  for(unsigned int num=0; b_it != b_to; ++b_it, ++num){
    // do temperature coupling for that bath?
    if (b_it->tau != -1 || immediate){

      DEBUG(7, "pre-scale ekin: " << b_it->ekin);

      double free_temp = 0.0;

      // small flexible constraints hack!
      if (sim.param().constraint.solute.algorithm == simulation::constr_flexshake){
	free_temp = 2 *
	  (b_it->ekin - conf.special().flexible_constraint.flexible_ekin[num])
	  / (b_it->dof * math::k_Boltzmann);
      }
      else{
	free_temp = 2 * 
	  b_it->ekin / (b_it->dof * math::k_Boltzmann);
      }
      
      // divide by zero measure...
      if (free_temp < math::epsilon) free_temp = b_it->temperature;

      // 2nd divide by zero measure... (if reference temperature is zero)
      if (free_temp < math::epsilon)
	b_it->scale = 1;
      else if (immediate)
	b_it->scale = sqrt(b_it->temperature / free_temp);
      else
	b_it->scale = sqrt(1.0 + sim.time_step_size() / b_it->tau *
			   (b_it->temperature / free_temp - 1));

      DEBUG(8, "free T " << free_temp << " dof " << b_it->dof);
      DEBUG(8, "ref. T " << b_it->temperature);
      DEBUG(8, "bath " << num << " scaling: " << b_it->scale);
      
    } // scaling ?
    else
      b_it->scale = 1;
  } // loop over the baths
  
}
