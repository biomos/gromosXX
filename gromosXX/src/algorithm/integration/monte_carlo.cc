/**
 * @file monte_carlo.cc
 * contains the implementation
 * for the Monte_Carlo class
 */

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <interaction/interaction.h>
#include <interaction/forcefield/forcefield.h>

#include "monte_carlo.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE integration

/**
 * Monte-Carlo step.
 */
int algorithm::Monte_Carlo::apply
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation &sim
 )
{
  if (sim.steps() == 0) return 0;
  
  if ((sim.steps() % sim.param().montecarlo.steps) == 0){
    move(topo, conf, sim);
    // reevaluate the forces
    m_ff.apply(topo, conf, sim);
    // see whether we accept
    accept(topo, conf, sim);
  }

  return 0;
}

// HARDCODE HARDCODE HARDCODE
const double q_max = 0.30;
const double q_min = -0.30;
const int atom = 0;
const double dq = 0.10;
const double T = 300;
// EDOCDRAH EDOCDRAH EDOCDRAH

double q_old;

/**
 * Monte-Carlo move.
 */
int algorithm::Monte_Carlo::move
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation &sim
 )
{
  conf.exchange_state();
  conf.current().box = conf.old().box;

  const int num_atoms = topo.num_atoms();
  for(int i=0; i < num_atoms; ++i){
    conf.current().vel(i) = conf.old().vel(i);
    conf.current().pos(i) = conf.old().pos(i);
    conf.current().force(i) = conf.old().force(i);
  }

  const int mc_step = sim.steps() / sim.param().montecarlo.steps;

  int sign = 1;
  if (false){
    sign = (mc_step % 2) ? 1 : -1;
  }
  else{
    double r =  gsl_rng_uniform(m_rng);
    if (r < 0.5) sign = -sign;
  }
  
  double q_new = topo.charge(atom) +  sign * dq;
  q_old = topo.charge(atom);

  // wrap around...
  if (q_new > q_max) q_new = q_min;
  if (q_new < q_min) q_new = q_max;
  
  std::cout << "MC: trying " << topo.charge(atom) << " -> " << q_new << std::endl;
  topo.charge()[atom] = q_new;

  return 0;
}

int algorithm::Monte_Carlo::accept
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation &sim
 )
{
  const int mc_step = (sim.steps()-1) / sim.param().montecarlo.steps;
  // const int sign = (mc_step % 2) ? 1 : -1;

  // get probability
  double delta = 0;
  const double beta = 1.0 / (math::k_Boltzmann * T);

  conf.current().energies.calculate_totals();
  conf.old().energies.calculate_totals();
  
  delta =
    beta * (conf.current().energies.potential_total - conf.old().energies.potential_total);

  std::cout << "MC: delta = " << delta << std::endl;

  double probability = 1.0;
  if (delta > 0.0)
    probability = exp(-delta);
  
  std::cout << "MC: old E_pot = " << conf.old().energies.potential_total
	    << "\tnew E_pot = " << conf.current().energies.potential_total << std::endl;
  std::cout << "MC: probability = " << probability << std::endl;
  
  const double r = gsl_rng_uniform(m_rng);
  if (r < probability){
    std::cout << "MC: switch succeeded!" << std::endl;
  }
  else{
    std::cout << "MC: reverting to charge " << q_old << std::endl;
    topo.charge()[atom] = q_old;
    // get old state back!
    conf.exchange_state();
  }

  std::cout << "MC: " << std::setw(7) << mc_step << "\tcharge " << std::setw(18) << topo.charge()[atom] << std::endl;
  return 0;
}
