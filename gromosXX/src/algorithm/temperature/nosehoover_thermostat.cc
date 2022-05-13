/**
 * @file nosehoover_thermostat.cc
 * methods of the Nose-Hoover Thermostat
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../configuration/state_properties.h"

#include "nosehoover_thermostat.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE

#define MODULE algorithm
#define SUBMODULE temperature

int algorithm::NoseHoover_Thermostat::init
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 std::ostream & os,
 bool quiet
 )
{
  if (sim.param().multibath.algorithm > 0){

    if (!quiet){
      if (sim.param().multibath.algorithm == 1){
	std::cout << "\tNose-Hoover temperature coupling\n";
      }
      else{
	std::cout << "\tNose-Hoover-Chain temperature coupling: using "
		  << sim.param().multibath.algorithm << " instances\n";
      }
    }
    
    std::vector<simulation::bath_struct>::iterator
      b_it = sim.multibath().begin(),
      b_to = sim.multibath().end();
  
    for( ; b_it != b_to; ++b_it){
      b_it->zeta.resize(sim.param().multibath.algorithm, 0.0);
    }
  }
  return 0;
}

int algorithm::NoseHoover_Thermostat::apply
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim
 )
{
  m_timer.start();

  assert(sim.param().multibath.algorithm > 0);
  
  if (sim.param().multibath.algorithm == 1)
    calc_scaling(topo, conf, sim);
  else
    calc_chain_scaling(topo, conf, sim);

  //--------------------------------
  // now we have the scaling factors
  //--------------------------------

  scale(topo, conf, sim);

  m_timer.stop();

  return 0;
  
}

void algorithm::NoseHoover_Thermostat
::calc_scaling
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim
 )
{
  DEBUG(7, "NoseHoover Temperature Scaling!!!!");

  // std::cout.precision(9);
  // std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
  
  // loop over the baths
  std::vector<simulation::bath_struct>::iterator
    b_it = sim.multibath().begin(),
    b_to = sim.multibath().end();
  
  for(unsigned int num=0; b_it != b_to; ++b_it, ++num){
    // do temperature coupling for that bath?
    if (b_it->tau != -1){

      DEBUG(7, "pre-scale ekin: " << b_it->ekin);

      double free_temp = 0.0;

      // small flexible constraints hack!
      if (sim.param().constraint.solute.algorithm == simulation::constr_flexshake){
	free_temp = 2 *
	  (b_it->ekin - conf.special().flexible_constraint.flexible_ekin[num]) / (b_it->dof * math::k_Boltzmann);
      }
      else{
	free_temp = 2 * 
	  b_it->ekin / (b_it->dof * math::k_Boltzmann);
      }
      
      // divide by zero measure...
      if (free_temp < math::epsilon) free_temp = b_it->temperature;

      b_it->zeta[0] += sim.time_step_size() / (b_it->tau * b_it->tau) * (free_temp / b_it->temperature - 1.0);
      b_it->scale = 1.0 - b_it->zeta[0] * sim.time_step_size();
      
      DEBUG(8, "free T " << free_temp << " dof " << b_it->dof);
      DEBUG(8, "ref. T " << b_it->temperature);
      DEBUG(8, "bath " << num << " scaling: " << b_it->scale
	    << " zeta: " << b_it->zeta[0]);
      
    } // scaling ?
    else
      b_it->scale = 1;
  } // loop over the baths
  
}

void algorithm::NoseHoover_Thermostat
::calc_chain_scaling
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim
 )
{
  DEBUG(7, "NoseHoover chains Temperature Scaling!!!!");

  // std::cout.precision(9);
  // std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
  
  const double dt = sim.time_step_size();

  // loop over the baths
  std::vector<simulation::bath_struct>::iterator
    b_it = sim.multibath().begin(),
    b_to = sim.multibath().end();
  
  for(unsigned int num=0; b_it != b_to; ++b_it, ++num){
    // do temperature coupling for that bath?
    if (b_it->tau != -1){

      DEBUG(7, "pre-scale ekin: " << b_it->ekin);

      double free_temp = 0.0;

      // small flexible constraints hack!
      if (sim.param().constraint.solute.algorithm == simulation::constr_flexshake){
	free_temp = 2 *
	  (b_it->ekin - conf.special().flexible_constraint.flexible_ekin[num]) /
	  (b_it->dof * math::k_Boltzmann);
      }
      else{
	free_temp = 2 * 
	  b_it->ekin / (b_it->dof * math::k_Boltzmann);
      }
      
      // divide by zero measure...
      if (free_temp < math::epsilon) free_temp = b_it->temperature;

      const int nhc = sim.param().multibath.algorithm;
      
      std::vector<double> tau(nhc, b_it->tau * b_it->tau / b_it->dof);
      tau[0] = b_it->tau * b_it->tau;

      assert(nhc > 1);
      
      b_it->zeta[nhc-1] += (tau[nhc-2] * b_it->zeta[nhc-2] * b_it->zeta[nhc-2]
                         - 1.0 / b_it->dof) / tau[nhc-1] * dt;
	
      for (int i = nhc - 2; i >= 1; i--){

	b_it->zeta[i] += ((tau[i-1] * b_it->zeta[i-1] * b_it->zeta[i-1] - 1.0 / b_it->dof)
			 / tau[i]  - b_it->zeta[i] * b_it->zeta[i+1]) * dt;
      }

      b_it->zeta[0] += ((free_temp / b_it->temperature - 1.0) / tau[0]
                         - b_it->zeta[1] * b_it->zeta[0] ) * dt;

      b_it->scale = 1.0 - b_it->zeta[0] * dt;
      
      DEBUG(8, "free T " << free_temp << " dof " << b_it->dof);
      DEBUG(8, "ref. T " << b_it->temperature);
      DEBUG(8, "bath " << num << " scaling: " << b_it->scale);
      
    } // scaling ?
    else
      b_it->scale = 1;
  } // loop over the baths
  
}
