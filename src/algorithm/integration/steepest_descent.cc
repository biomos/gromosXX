/**
 * @file steepest_descent.cc
 * contains the implementation
 * for steepest descent energy minimisation
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "steepest_descent.h"

#include "../../util/error.h"

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
       simulation::Simulation &sim,
       std::ostream & os,
       bool quiet)
{
  if (!quiet) {
    os << "ENERGY MINIMISATION\n"
       << "\tsteepest descent\n"
       << "\tminimum energy criterion : " << sim.param().minimise.dele << "\n"
       << "\tstarting step size       : " << sim.param().minimise.dx0 << "\n"
       << "\tmaximum step size        : " << sim.param().minimise.dxm << "\n"
       << "\tminimum steps            : " << sim.param().minimise.nmin << "\n";
  }
  
  if (sim.param().minimise.flim != 0) {
    if (!quiet) {
      os << "\tlimiting the force to    : "
	 << sim.param().minimise.flim << "\n";
    }
  }
  
  if (!quiet) {
    os << "END\n";
  }
  
  conf.old().vel = 0.0;
  conf.current().vel = 0.0;
  
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
  if (sim.steps() > unsigned(sim.param().minimise.nmin)) {
    // check whether minimum reached...
    conf.current().energies.calculate_totals();

    const double eold = conf.old().energies.potential_total + conf.old().energies.special_total;
    const double ecur = conf.current().energies.potential_total +  conf.current().energies.special_total;

    DEBUG(10, "epot + espec: " 
	  << ecur
	  << "\told epot + espec: " 
	  << eold
	  << "\tdiff: " << fabs(ecur - eold)
	  << "\tdele: " << sim.param().minimise.dele);
    
    if (fabs(ecur - eold) < sim.param().minimise.dele) {
      std::cout << "STEEPEST DESCENT:\tMINIMUM REACHED\n";
      return E_MINIMUM_REACHED;
    }
    
    if (ecur < eold) {
      sim.minimisation_step_size() *= 1.2;
      if (sim.minimisation_step_size() > sim.param().minimise.dxm) {
        sim.minimisation_step_size() = sim.param().minimise.dxm;
      }
    }
    else {
      sim.minimisation_step_size() *= 0.5;
    }
  }
  else {
    sim.minimisation_step_size() = sim.param().minimise.dx0;
  }

  // limit the maximum force!
  if (sim.param().minimise.flim != 0.0) {
    for(unsigned int i=0; i<topo.num_atoms(); ++i) {
      const double fs = math::abs(conf.current().force(i));
      if (fs > sim.param().minimise.flim) {
        DEBUG(15,"fs = " << fs 
             << ", flim = " << sim.param().minimise.flim);
        DEBUG(15,"force (" << i << ") = " << math::v2s(conf.current().force(i))
             << ", factor = " << sim.param().minimise.flim / fs);
        conf.current().force(i) *= sim.param().minimise.flim / fs;
      }
    }
  }

  // <f|f>^-0.5
  // double f = math::sum(math::abs2(conf.current().force));
  double f = 0.0;
  for(unsigned int i=0; i<topo.num_atoms(); ++i) {
    f += math::abs2(conf.current().force(i));
  }
  
  if (f < math::epsilon) {
    f = 1.0;
  }
  else {
    f = 1.0 / sqrt(f);
  }

/*
#ifdef HAVE_ISNAN
  if (std::isnan(f)){
    io::messages.add("force is NaN", "Steepest_Descent", io::message::error);
    return E_NAN;
  }
#endif
*/
  if(math::gisnan(f)) {
    io::messages.add("force is NaN", "Steepest_Descent", io::message::error);
    return E_NAN;
  }
  
  conf.exchange_state();

  for(unsigned int i=0; i<topo.num_atoms(); ++i) {
    conf.current().pos(i) = conf.old().pos(i) + sim.minimisation_step_size() * f * conf.old().force(i);
  }

  conf.old().vel = 0.0;
  conf.current().vel = 0.0;
  
  return 0;
}

