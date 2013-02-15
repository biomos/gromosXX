/**
 * @file ramd_interaction.cc
 * template methods of ramd_interaction
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"
#include "../../interaction/forcefield/forcefield.h"

#include "../../math/periodicity.h"

// special interactions
#include "../../interaction/interaction_types.h"

#include "../../interaction/special/ramd_interaction.h"

#include "create_special.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

/**
 * calculate_interactions
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_ramd_interactions
(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 gsl_rng * m_rng)
{

  // check for update force-direction
  if(!(sim.steps() % sim.param().ramd.steps))
  {
    DEBUG(7, "Considering new force at step " << sim.steps());
  
    math::Periodicity<B> periodicity(conf.current().box);
    // check for sufficient movement
    // calculate com (is everything gathered?)
    math::Vec r;
    math::Vec com(0.0,0.0,0.0);

    // calculate some properties for time averaging
    double eterm=0.0;
    double mem_dec=1.0;
    if(sim.param().ramd.do_ta){
      eterm = exp( -sim.param().ramd.steps*sim.param().step.dt / 
		   sim.param().ramd.tau);
      mem_dec = 1.0 - eterm;
    }
    
    std::set<unsigned int>::const_iterator
      it = sim.param().ramd.atom.begin(),
      i0 = sim.param().ramd.atom.begin(),
      to = sim.param().ramd.atom.end();

    for(; it!=to; ++it){
      periodicity.nearest_image(conf.current().pos(*it), 
				conf.current().pos(*i0), r);
      com += topo.mass()(*it) * r;
    }
    
    com /= conf.special().ramd.total_mass;
    com += conf.current().pos(*i0);
    
    double r_min2 = sim.param().ramd.r_min * sim.param().ramd.r_min;
    periodicity.nearest_image(com, conf.special().ramd.old_com, r);
    double r2 = abs2(r);
    double dist = sqrt(r2);
    
    if(sim.param().ramd.do_ta && sim.steps()!=0 )
      conf.special().ramd.ta_average = mem_dec * dist + 
	eterm * conf.special().ramd.ta_average;
    

    if(r2 < r_min2){
      
      // obtain new random direction
      r(0) = gsl_rng_uniform(m_rng)-0.5;
      r(1) = gsl_rng_uniform(m_rng)-0.5;
      r(2) = gsl_rng_uniform(m_rng)-0.5;
      
      conf.special().ramd.force_direction = r / abs(r);

      DEBUG(7,"new force direction " 
	    << conf.special().ramd.force_direction(0) << " " 
	    << conf.special().ramd.force_direction(1) << " " 
	    << conf.special().ramd.force_direction(2));
      
    }
    if(conf.special().ramd.ta_average < sim.param().ramd.ta_min){
      topo.lambda(1.0);
    }
    else{
      topo.lambda(0.0);
    }
    DEBUG(7, "lambda value after RAMD " << topo.lambda());
    
    topo.update_for_lambda();
   
    // store com
    conf.special().ramd.old_com = com;
    
  }
     
  // apply random force
  std::set<unsigned int>::const_iterator
    it = sim.param().ramd.atom.begin(),
    to = sim.param().ramd.atom.end();
  
  for(; it!=to; ++it)
    conf.current().force(*it) +=
      sim.param().ramd.fc * conf.special().ramd.force_direction * 
      topo.mass()(*it)/conf.special().ramd.total_mass;
    
  return 0;
}

/**
 * init
 */
int interaction::RAMD_Interaction::init
(
 topology::Topology &topo, 
 configuration::Configuration &conf,
 simulation::Simulation &sim,
 std::ostream &os,
 bool quiet) 
{
  if (!quiet){
    os << "\nRAMD: Random Acceleration MD simulation\n";
    os << "\tForce constant: " << sim.param().ramd.fc << "\n"
       << "\tNumber of steps between updates: " << sim.param().ramd.steps << "\n"
       << "\tMinimum distance for COM to travel: " << sim.param().ramd.r_min << "\n"
       << "\tWrite COM every " << sim.param().ramd.every << " steps\n"
       << "\tRAMD atoms:\n\t";
    std::set<unsigned int>::iterator
      it = sim.param().ramd.atom.begin(),
      to = sim.param().ramd.atom.end();
    for(int i=1;it!=to; ++it, ++i){
      os << std::setw(6) << *it+1;
      if(i%10==0) os << "\n\t";
    }
  }
  
  // calculate total mass of ramd atoms
  conf.special().ramd.total_mass = 0.0;
  
  std::set<unsigned int>::const_iterator
    it = sim.param().ramd.atom.begin(),
    to = sim.param().ramd.atom.end();
  for( ;it!=to; ++it)
    conf.special().ramd.total_mass += topo.mass()(*it);
  if(!quiet)
    os << "\n\n\tTotal mass of RAMD atoms: " << conf.special().ramd.total_mass
       << "\n\nEND\n";

  // come up with initial direction
  math::Vec r;
  // initialize random number seed
  gsl_rng_set(m_rng, sim.param().start.ig);
  r(0) = gsl_rng_uniform(m_rng) - 0.5;
  r(1) = gsl_rng_uniform(m_rng) - 0.5;
  r(2) = gsl_rng_uniform(m_rng) - 0.5;
  
  conf.special().ramd.force_direction = r / abs(r);

  
  return 0;
};



int interaction::RAMD_Interaction::calculate_interactions
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim)
{
  SPLIT_VIRIAL_BOUNDARY(_calculate_ramd_interactions,
                        topo, conf, sim, m_rng);

  return 0;
}

