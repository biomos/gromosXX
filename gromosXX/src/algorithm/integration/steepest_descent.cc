/**
 * @file steepest_descent.cc
 * contains the implementation
 * for steepest descent energy minimisation
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include "steepest_descent.h"

#include <util/error.h>

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE integration

/**
 * steepest descent step.
 */
int algorithm::Steepest_Descent
::init(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation &sim)
{
  std::cout << "ENERGY MINIMISATION\n"
	    << "\tsteepest descent\n"
	    << "\tminimum energy criterion : " << sim.param().minimise.dele << "\n"
	    << "\tstarting step size       : " << sim.param().minimise.dx0 << "\n"
	    << "\tmaximum step size        : " << sim.param().minimise.dxm << "\n"
	    << "\tminimum steps            : " << sim.param().minimise.nmin << "\n";
  
  if (sim.param().minimise.flim != 0)
    std::cout << "\tlimiting the force to    : "
	      << sim.param().minimise.flim << "\n";
  
  std::cout << "END\n";
  return 0;
}

/**
 * steepest descent step.
 */
int algorithm::Steepest_Descent
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation &sim)
{
  // only if it's not the first step
  if (sim.steps() > unsigned(sim.param().minimise.nmin)){
    // check whether minimum reached...
    conf.current().energies.calculate_totals();

    DEBUG(10, "epot: " << conf.current().energies.potential_total
	  << "\told epot: " << conf.old().energies.potential_total
	  << "\tdiff: " << fabs(conf.current().energies.potential_total -
				conf.old().energies.potential_total)
	  << "\tdele: " << sim.param().minimise.dele);
    
    if (fabs(conf.current().energies.potential_total -
	     conf.old().energies.potential_total) < sim.param().minimise.dele)
    {
      std::cout << "STEEPEST DESCENT:\tMINIMUM REACHED\n";
      return E_MINIMUM_REACHED;
    }
    
    if (conf.current().energies.potential_total < conf.old().energies.potential_total){
      sim.minimisation_step_size() *= 1.2;
      if (sim.minimisation_step_size() > sim.param().minimise.dxm)
	sim.minimisation_step_size() = sim.param().minimise.dxm;
    }
    else
      sim.minimisation_step_size() *= 0.5;
  }
  else
    sim.minimisation_step_size() = sim.param().minimise.dx0;

  // limit the maximum force!
  if (sim.param().minimise.flim != 0.0){
    for(unsigned int i=0; i<topo.num_atoms(); ++i){
      const double fs = math::dot(conf.current().force(i), 
				  conf.current().force(i));
      if (fs > sim.param().minimise.flim)
	conf.current().force(i) *= sim.param().minimise.flim / fs;
    }
  }

  // <f|f>^-0.5
  double f = math::sum(math::dot(conf.current().force, conf.current().force));
  f = 1.0 / sqrt(f);

#ifdef HAVE_ISNAN
  if (isnan(f)){
    io::messages.add("force is NaN", "Steepest_Descent", io::message::error);
    return E_NAN;
  }
#endif

  conf.exchange_state();

  conf.current().pos = conf.old().pos + sim.minimisation_step_size() * f *
    conf.old().force;

  conf.old().vel = 0.0;
  conf.current().vel = 0.0;
  
  return 0;
}

