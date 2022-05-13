/**
 * @file eds_distance_restraint_interaction.cc
 * template methods of Eds_Distance_Restraint_Interaction
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"

#include "../../math/periodicity.h"

// special interactions
#include "../../interaction/interaction_types.h"

#include "../../interaction/special/eds_distance_restraint_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

/**
 * calculate position restraint interactions
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_eds_distance_restraint_interactions
(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim, double exponential_term)
{
  // loop over the position restraints
  std::vector<topology::eds_distance_restraint_struct>::const_iterator 
    it = topo.eds_distance_restraints().begin(),
    to = topo.eds_distance_restraints().end();
  
  
  // math::VArray &pos   = conf.current().pos;
  // math::VArray &force = conf.current().force;
  math::Vec v, f;

  double energy = 0.0;

  math::Periodicity<B> periodicity(conf.current().box);

  // for(int i=0; it != to; ++it, ++i){
  for(; it != to; ++it){

    periodicity.nearest_image(it->v1.pos(conf,topo), it->v2.pos(conf,topo), v);
#ifndef NDEBUG
    for(int i=0;i<it->v1.size();i++){
      DEBUG(10, "v1 (atom " << i+1 << "): " << it->v1.atom(i)+1);
    }
    for(int i=0;i<it->v2.size();i++){
      DEBUG(10, "v2 (atom " << i+1 << "): " << it->v2.atom(i)+1);
    }
#endif
    DEBUG(10, "pos(v1) = " << math::v2s(it->v1.pos(conf,topo)));
    DEBUG(10, "pos(v2) = " << math::v2s(it->v2.pos(conf,topo)));

    DEBUG(9, "DISTANCERES v : " << math::v2s(v));
    
    double dist = math::abs(v);
    const unsigned int numstates = sim.param().eds.numstates;
    // loop over all eds states
    for(unsigned int state = 0; state < numstates; state++){
      DEBUG(9, "DISTANCERES dist : " << dist << " r0 [" << state << "]" << it->r0[state]);
    
      if(it->rah*dist < it->rah * it->r0[state]){
        DEBUG(9, "DISTANCERES  : restraint fulfilled");
        f=0;
      }
      else if(fabs(it->r0[state] - dist) < sim.param().distanceres.r_linear){
        DEBUG(9, "DISTANCERES  : harmonic");
        f = - sim.param().distanceres.K * (dist - it->r0[state]) * v / dist;
      }
      else{
        DEBUG(9, "DISTANCERES  : linear");
        if(dist < it->r0[state])
          f =  sim.param().distanceres.r_linear * sim.param().distanceres.K  * v / dist;
        else
          f = -sim.param().distanceres.r_linear * sim.param().distanceres.K  * v / dist;
      }
      
      if(abs(sim.param().distanceres.distanceres) == 1)
        ;
      else if(abs(sim.param().distanceres.distanceres) == 2)
        f *= it->w0[state];
      else
        f = 0;
      
      DEBUG(9, "Distanceres force : " << math::v2s(f));
      
      it->v1.force(conf, topo,  f, conf.special().eds.force_endstates[state]);
      it->v2.force(conf, topo, -f, conf.special().eds.force_endstates[state]);
      
      if(it->rah * dist < it->rah * it->r0[state])
        energy = 0;
      else if(fabs(it->r0[state] - dist) < sim.param().distanceres.r_linear)
        energy = 0.5 * sim.param().distanceres.K * (dist - it->r0[state]) * (dist - it->r0[state]);
      else{
        if(dist < it->r0[state])
          energy = -sim.param().distanceres.K * sim.param().distanceres.r_linear *
                  (dist - it->r0[state] + 0.5 * sim.param().distanceres.r_linear);
        else
          energy =  sim.param().distanceres.K * sim.param().distanceres.r_linear *
                  (dist - it->r0[state] - 0.5 * sim.param().distanceres.r_linear);
      }
      
      if(abs(sim.param().distanceres.distanceres) == 1)
        ;
      else if(abs(sim.param().distanceres.distanceres) == 2)
        energy *= it->w0[state];
      else
        energy=0;
      
      DEBUG(9, "Distanceres energy : " << energy);
      conf.current().energies.eds_vi[state] += energy;
      conf.current().energies.eds_vi_special[state] += energy;
      
      // std::cout.precision(9);
      // std::cout.setf(std::ios::fixed, std::ios::floatfield);
      // std::cout << "DISTANCERES_" << i << ": "
      // << dist << "\t\t" << math::dot(f, v) / dist
      // << "\t\t" << energy << "\n";
    } // loop over eds states
  }
  
  return 0;
}

int interaction::Eds_Distance_Restraint_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  SPLIT_VIRIAL_BOUNDARY(_calculate_eds_distance_restraint_interactions,
			topo, conf, sim, exponential_term);
  
  return 0;
}

int interaction::Eds_Distance_Restraint_Interaction::init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os,
		     bool quiet) 
{
 
  if (!quiet) {
    os << "eds perturbed distance restraint interaction";
    os << std::endl;
  }
  return 0;
}
