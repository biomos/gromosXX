/**
 * @file monte_carlo.cc
 * contains the implementation
 * for the Monte_Carlo class
 */

#ifdef XXMPI
#include <mpi.h>
#endif

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../interaction/interaction.h"
#include "../../interaction/forcefield/forcefield.h"

#include "monte_carlo.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE integration

/**
 * Chemical Monte-Carlo step.
 */
int algorithm::Monte_Carlo::apply(
topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation &sim
) {
  DEBUG(12, "Chemical Monte-Carlo apply");
  if (sim.steps() == 0) return 0;
  // do a CMC step?
  if ((sim.steps() % sim.param().montecarlo.steps) == 0){
    int rank=0;
#ifdef XXMPI
    if (sim.mpi) {
      rank = MPI::COMM_WORLD.Get_rank();
    }
#endif
    DEBUG(14,"Chemical Monte-Carlo apply, rank = " << rank);
    // MASTER
    if(rank == 0){
      move(topo, conf, sim);
      // reevaluate the forces
      DEBUG(14, "Chemical Monte-Carlo apply, MASTER");
      m_ff.apply(topo, conf, sim);
    }
    
#ifdef XXMPI
    // SLAVE
    
    else{
      DEBUG(14, "Chemical Monte-Carlo apply, SLAVE");
      interaction::Interaction * nb = m_ff.interaction("NonBonded");
      if (nb == NULL){
        std::cerr << "MPI slave: could not get NonBonded interactions from forcefield"
        << "\n\t(internal error)"
        << std::endl;
        MPI::Finalize();
        return 1;
      }
 
      if ((nb->calculate_interactions(topo, conf, sim)) != 0){
        std::cerr << "MPI slave " << rank << ": error in nonbonded calculation!\n" << std::endl;
      }
    }
    
#endif
    if(rank==0){
      // see whether we accept
      accept(topo, conf, sim);
    }

  }

  return 0;
}

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
  DEBUG(12,"Chemical Monte-Carlo move");
  conf.exchange_state();
  conf.current().box = conf.old().box;

  const int num_atoms = topo.num_atoms();
  for(int i=0; i < num_atoms; ++i){
    conf.current().vel(i) = conf.old().vel(i);
    conf.current().pos(i) = conf.old().pos(i);
    conf.current().force(i) = conf.old().force(i);
  }
  
  // This function returns a double precision floating
  // point number uniformly distributed in the range [0,1). 
  // The range includes 0.0 but excludes 1.0.
  // should we have [0,1]?
  double ran = gsl_rng_uniform(m_rng)-0.5;
  double lambda_new = topo.lambda() + ran*sim.param().montecarlo.dlambda;
  m_lambda_old=topo.lambda();
 
  // wrap around
  if(lambda_new < 0){
     lambda_new = lambda_new - static_cast<int>(lambda_new) + 1;
  }
  else{
    lambda_new = lambda_new - static_cast<int>(lambda_new);
  }
  DEBUG(14,"CMC: trying lambda " << topo.lambda() << " -> " << lambda_new);
  // change the lambda value
  topo.lambda(lambda_new);
  // update masses
  topo.update_for_lambda();
  
  return 0;
}

int algorithm::Monte_Carlo::accept
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation &sim
 )
{
  DEBUG(12,"Chemcial Monte-Carlo accept");
  const int mc_step = (sim.steps()-1) / sim.param().montecarlo.steps;
    
  // get probability
  double delta = 0;
  // change T to the actual temperature!!!
  // --> make sure we only have one T! (in_parameter.cc)
  //sim.param().multibath.multibath.bath(0)
  double temperature = sim.param().multibath.multibath.bath(0).temperature;
  const double beta = 1.0 / (math::k_Boltzmann * temperature);
  
  conf.current().energies.calculate_totals();
  conf.old().energies.calculate_totals();
  
  double probability = 1.0;
  
  delta =beta *
          (conf.current().energies.potential_total
          - conf.old().energies.potential_total);
  
   DEBUG(12, "CMC: beta*(Epot(new)-Epot(old)) = " << delta );
  
  
  /*
   * if Epot(new) > Epot(old) only accept with probability exp(-delta)
   * else if Epot(new) < Epot(old) accept always (prob.=1.0)
   */
  
  if (delta > 0.0) probability = exp(-delta);
  
  DEBUG(12, "CMC: old E_pot = " << conf.old().energies.potential_total
           << "\tnew E_pot = " << conf.current().energies.potential_total);
  DEBUG(12, "CMC: probability = " << probability);
 
  const double r = gsl_rng_uniform(m_rng);
  if (r < probability){
    DEBUG(12, "CMC: switch succeeded!");
  }
  else{
    DEBUG(12,"CMC: reverting to lambda " << m_lambda_old);
    // change the lambda value back
    topo.lambda(m_lambda_old);
    // update masses
    topo.update_for_lambda();
    
    // get old state back! - does this work?
    conf.exchange_state();
  }
  DEBUG(12, "MC: " << std::setw(7) << mc_step << "\tlambda " << std::setw(18) <<topo.lambda() );
  return 0;
}
