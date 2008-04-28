/**
 * @file distance_restraint_interaction.cc
 * template methods of Distance_Restraint_Interaction
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>
#include <interaction/interaction.h>

#include <math/periodicity.h>

// special interactions
#include <interaction/interaction_types.h>

#include <interaction/special/distance_restraint_interaction.h>

#include <util/template_split.h>
#include <util/debug.h>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

/**
 * calculate position restraint interactions
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_distance_restraint_interactions
(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim, double exponential_term)
{
  // loop over the position restraints
  std::vector<topology::distance_restraint_struct>::const_iterator 
    it = topo.distance_restraints().begin(),
    to = topo.distance_restraints().end();
  
  std::vector<double>::iterator ave_it = conf.special().distanceres_av.begin();

  // math::VArray &pos   = conf.current().pos;
  // math::VArray &force = conf.current().force;
  math::Vec v, f;

  double energy;

  math::Periodicity<B> periodicity(conf.current().box);

  // for(int i=0; it != to; ++it, ++i){
  for(; it != to; ++it, ++ave_it){

    periodicity.nearest_image(it->v1.pos(conf), it->v2.pos(conf), v);
#ifndef NDEBUG
    for(int i=0;i<it->v1.size();i++){
      DEBUG(10, "v1 (atom " << i+1 << "): " << it->v1.atom(i)+1);
    }
    for(int i=0;i<it->v2.size();i++){
      DEBUG(10, "v2 (atom " << i+1 << "): " << it->v2.atom(i)+1);
    }
#endif
    DEBUG(10, "pos(v1) = " << math::v2s(it->v1.pos(conf)));
    DEBUG(10, "pos(v2) = " << math::v2s(it->v2.pos(conf)));

    DEBUG(9, "DISTANCERES v : " << math::v2s(v));
    
    double dist = math::abs(v);
    DEBUG(9, "DISTANCERES dist : " << dist << " r0 : " << it->r0);
    if (sim.param().distanceres.distanceres < 0) {
      (*ave_it) = (1.0 - exponential_term) * pow(dist, -3.0) + 
                   exponential_term * (*ave_it);
      dist = pow(*ave_it, -1.0 / 3.0);
      DEBUG(9, "DISTANCERES average dist : " << dist);
    }
      

    if(it->rah*dist < it->rah * it->r0){
      DEBUG(9, "DISTANCERES  : restraint fulfilled");
      f=0;
    }
    else if(fabs(it->r0 - dist) < sim.param().distanceres.r_linear){
      DEBUG(9, "DISTANCERES  : harmonic");
      f = - sim.param().distanceres.K * (dist - it->r0) * v / dist;
    }
    else{
      DEBUG(9, "DISTANCERES  : linear");
      if(dist < it->r0)
	f =  sim.param().distanceres.r_linear * sim.param().distanceres.K  * v / dist;
      else
	f = -sim.param().distanceres.r_linear * sim.param().distanceres.K  * v / dist;
    }
    
    if(abs(sim.param().distanceres.distanceres) == 1)
      ;
    else if(abs(sim.param().distanceres.distanceres) == 2)
      f *= it->w0;
    else 
      f = 0;
    
    DEBUG(9, "Distanceres force : " << math::v2s(f));

    it->v1.force(conf,  f);
    it->v2.force(conf, -f);  

    if(it->rah * dist < it->rah * it->r0)
      energy = 0;
    else if(fabs(it->r0 - dist) < sim.param().distanceres.r_linear)
      energy = 0.5 * sim.param().distanceres.K * (dist - it->r0) * (dist - it->r0);
    else{
      if(dist < it->r0)
	energy = -sim.param().distanceres.K * sim.param().distanceres.r_linear *
	  (dist - it->r0 + 0.5 * sim.param().distanceres.r_linear);
      else
	energy =  sim.param().distanceres.K * sim.param().distanceres.r_linear *
	  (dist - it->r0 - 0.5 * sim.param().distanceres.r_linear);
    }
    
    if(abs(sim.param().distanceres.distanceres) == 1)
      ;
    else if(abs(sim.param().distanceres.distanceres) == 2)
     energy *= it->w0;
    else
      energy=0;

    DEBUG(9, "Distanceres energy : " << energy);

    conf.current().energies.distanceres_energy[topo.atom_energy_group()
					    [it->v1.atom(0)]] += energy;
    
    // std::cout.precision(9);
    // std::cout.setf(std::ios::fixed, std::ios::floatfield);
    // std::cout << "DISTANCERES_" << i << ": " 
    // << dist << "\t\t" << math::dot(f, v) / dist 
    // << "\t\t" << energy << "\n";
    
  }
  
  return 0;
}

int interaction::Distance_Restraint_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  SPLIT_VIRIAL_BOUNDARY(_calculate_distance_restraint_interactions,
			topo, conf, sim, exponential_term);
  
  return 0;
}

/**
 * calculate position restraint interactions
 */
template<math::boundary_enum B>
static void _init_averages
(topology::Topology & topo,
 configuration::Configuration & conf)
{
  math::Periodicity<B> periodicity(conf.current().box);
  math::Vec v;
  
  for(std::vector<topology::distance_restraint_struct>::const_iterator
        it = topo.distance_restraints().begin(),
        to = topo.distance_restraints().end(); it != to; ++it) {
    periodicity.nearest_image(it->v1.pos(conf), it->v2.pos(conf), v);
    conf.special().distanceres_av.push_back(pow(math::abs(v), -3.0));
  }
}

int interaction::Distance_Restraint_Interaction::init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os,
		     bool quiet) 
{
  if (sim.param().distanceres.distanceres < 0) {
    exponential_term = std::exp(- sim.time_step_size() / 
                                  sim.param().distanceres.tau);
    
    if (!sim.param().distanceres.read) {
      // reset averages to r_0
      SPLIT_BOUNDARY(_init_averages, topo, conf);
    }
  }
  
  if (!quiet) {
    os << "Distance restraint interaction";
    if (sim.param().distanceres.distanceres < 0)
      os << "with time-averaging";
    os << std::endl;
  }
  return 0;
}
