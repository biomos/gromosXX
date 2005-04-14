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
 simulation::Simulation & sim)
{
  // loop over the position restraints
  std::vector<topology::distance_restraint_struct>::const_iterator 
    it = topo.distance_restraints().begin(),
    to = topo.distance_restraints().end();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec v, f;

  double energy;

  math::Periodicity<B> periodicity(conf.current().box);

  for( ; it != to; ++it){

    periodicity.nearest_image(it->v1.pos(conf), it->v2.pos(conf), v);

    DEBUG(10, "v1 : " << it->v1.atom(0) << " - " << it->v1.atom(1)
	  <<  " - " <<it->v1.atom(2) << " - " << it->v1.atom(3));
    DEBUG(10, "v2 : " << it->v2.atom(0) << " - " << it->v2.atom(1)
	  <<  " - " << it->v2.atom(2) << " - " << it->v2.atom(3));

    DEBUG(10, "pos(v1) = " << math::v2s(it->v1.pos(conf)));
    DEBUG(10, "pos(v2) = " << math::v2s(it->v2.pos(conf)));

    DEBUG(9, "DISTREST v : " << math::v2s(v));
    
    double dist = math::abs(v);

    DEBUG(9, "DISTREST dist : " << dist << " r0 : " << it->r0);

    if(it->rah*dist < it->rah * it->r0){
      DEBUG(9, "DISTREST  : restraint fulfilled");
      f=0;
    }
    else if(fabs(it->r0 - dist) < sim.param().distrest.r_linear){
      DEBUG(9, "DISTREST  : harmonic");
      f = - sim.param().distrest.K * (dist - it->r0) * v / dist;
    }
    else{
      DEBUG(9, "DISTREST  : linear");
      if(dist < it->r0)
	f =  sim.param().distrest.r_linear * sim.param().distrest.K  * v / dist;
      else
	f = -sim.param().distrest.r_linear * sim.param().distrest.K  * v / dist;
    }
    
    if(sim.param().distrest.distrest == 1)
      ;
    else if(sim.param().distrest.distrest == 2)
      f *= it->w0;
    else 
      f = 0;
    
    DEBUG(9, "Distrest force : " << math::v2s(f));

    it->v1.force(conf,  f);
    it->v2.force(conf, -f);  

    if(it->rah * dist < it->rah * it->r0)
      energy = 0;
    else if(fabs(it->r0 - dist) < sim.param().distrest.r_linear)
      energy = 0.5 * sim.param().distrest.K * (dist - it->r0) * (dist - it->r0);
    else{
      if(dist < it->r0)
	energy = -sim.param().distrest.K * sim.param().distrest.r_linear *
	  (dist - it->r0 + 0.5 * sim.param().distrest.r_linear);
      else
	energy =  sim.param().distrest.K * sim.param().distrest.r_linear *
	  (dist - it->r0 - 0.5 * sim.param().distrest.r_linear);
    }
    
    if(sim.param().distrest.distrest == 1)
      ;
    else if(sim.param().distrest.distrest == 2)
     energy *= it->w0;
    else
      energy=0;

    DEBUG(9, "Distrest energy : " << energy);

    conf.current().energies.distrest_energy[topo.atom_energy_group()
					    [it->v1.atom(0)]] += energy;
    
    std::cout.precision(9);
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout << "DISTREST " << dist << "\t\t" << abs(f) << "\t\t" << energy << "\n";
    
  }
  
  return 0;
}

int interaction::Distance_Restraint_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{

  SPLIT_VIRIAL_BOUNDARY(_calculate_distance_restraint_interactions,
			topo, conf, sim);
  
  return 0;
}
