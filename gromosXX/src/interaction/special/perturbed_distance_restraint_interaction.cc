/**
 * @file perturbed_distance_restraint_interaction.cc
 * template methods of Perturbed_Distance_Restraint_Interaction
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

#include <interaction/special/perturbed_distance_restraint_interaction.h>

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
static int _calculate_perturbed_distance_restraint_interactions
(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim)
{
  // loop over the perturbed distance restraints
  std::vector<topology::perturbed_distance_restraint_struct>::const_iterator 
    it = topo.perturbed_distance_restraints().begin(),
    to = topo.perturbed_distance_restraints().end();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec v, f;

  double energy, energy_derivativ;


  math::Periodicity<B> periodicity(conf.current().box);

  DEBUG(7, "perturbed distance restraint interactions : " 
	<< topo.perturbed_distance_restraints().size());

  for( ; it != to; ++it){

    periodicity.nearest_image(it->v1.pos(pos), it->v2.pos(pos), v);
 
    DEBUG(10, "v1 : " << it->v1.atom(0) << " - " << it->v1.atom(1)
	  <<  " - " <<it->v1.atom(2) << " - " << it->v1.atom(3));
    DEBUG(10, "v2 : " << it->v2.atom(0) << " - " << it->v2.atom(1)
	  <<  " - " << it->v2.atom(2) << " - " << it->v2.atom(3));
    
    DEBUG(9, "PERTDISTREST v : " << math::v2s(v));
   
    double dist = abs(v);

    const double l = topo.lambda();
    const double l2 = l*l;
    const double w0 = (1-l)*it->A_w0 + l*it->B_w0;
    const double r0 = (1-l)*it->A_r0 + l*it->B_r0;
    const double K = sim.param().distrest.K;
    const double r_l = sim.param().distrest.r_linear; 
    const double D_r0 = it->B_r0 - it->A_r0;
    DEBUG(9, "PERTDISTREST dist : " << dist << " r0 " << r0 << " rah " << it->rah);  

    if(it->rah*dist < it->rah*(r0))
      {
	DEBUG(9, "PERTDISTREST  : (rep / attr) restraint fulfilled");
	f=0*v;		
      }    
    else if(fabs(r0 - dist) < r_l)
      {
	DEBUG(9, "PERTDISTREST  :  harmonic");
	f = -4*( l*(1-l)*K*(dist - (r0))*v / dist); 
      }    
    else{
      DEBUG(9, "PERTDISTREST  : (rep / attr) linear");
      if(dist<r0)
	f =4*( l*(1-l) * r_l * K * v / dist);
      else
	f =-4*( l*(1-l) * r_l * K * v / dist);
    }
  
    if(sim.param().distrest.distrest == 1)
      ;      
    else if(sim.param().distrest.distrest == 2)
      f=f*w0;
    else
      f=f*0;

       
    it->v1.force(pos,f, force);
    it->v2.force(pos,-f, force);  
  
    if(it->rah*dist < it->rah*r0)
      {
	DEBUG(10, "PERTDISTREST: (rep / att ) restraint fulfilled");
      	energy = 0;
      }
    else if(fabs(r0 - dist) < r_l)
      {
	DEBUG(10, "PERTDISTREST: harmonic");	
	energy = 2 * l * (1-l) * K * (dist - r0)* (dist - r0);
      }
    else 
      {
	DEBUG(10, "PERTDISTREST: (rep / att ) linear");	
	if(dist<r0)
	  energy= -4*l*(1-l)* K *r_l *(dist  + 0.5 * r_l - r0);
	else
	  energy=  4*l*(1-l)* K *r_l *(dist  - 0.5 * r_l - r0);
      }
    
    
    if(sim.param().distrest.distrest == 1)
      ;
    else if(sim.param().distrest.distrest == 2)
      energy = energy*w0;
    else
      energy=0;
    
    
    DEBUG(10, "PERTDISTREST Energy : " << energy);
    conf.current().energies.distrest_energy[topo.atom_energy_group()
					    [it->v1.atom(0)]] += energy;
    
    if(it->rah*dist <it->rah*(r0))
      energy_derivativ = 0;
    
    else if(sim.param().distrest.distrest == 1)
      {
	if(fabs(r0-dist)<r_l)
	  energy_derivativ = 2*((1-2*l)*K*(dist - (r0)) * (dist - (r0)))
	    -4*(l*(1-l)*K*(dist - (r0))* (D_r0));
       	else 
	  {
	    if(dist<r0)
	      energy_derivativ= -4*r_l*K*((1-2*l)*(dist  +0.5 *r_l - (r0)) -l*(1-l)*(D_r0));
	    else
	      energy_derivativ= 4*r_l*K*((1-2*l)*(dist  - 0.5 *r_l - (r0)) -l*(1-l)*(D_r0));
	    }	
      }
    
    else if(sim.param().distrest.distrest == 2)
      {
	if(fabs(r0-dist)<r_l)
	  energy_derivativ = 2*(K*((1-4*l+3*l2)*it->A_w0 +(2*l-3*l2)*it->B_w0)*
				(dist - (r0))*(dist - (r0)))
	    -4*l*K*(((1-2*l+l2)*it->A_w0+ l*(1-l)*it->B_w0)*(dist - (r0)) * (D_r0));
	
	else 
	  {
	    if(dist<r0)
	      energy_derivativ=- 4*r_l*K*(((1-4*l+3*l2)*it->A_w0 + (2*l-3*l2)*it->B_w0)*(dist  + 0.5*r_l - (r0))
				     -l*((1-2*l+l2)*it->A_w0 +l*(1-l)*it->B_w0)*(D_r0));
	    else   
	      energy_derivativ= 4*r_l*K*(((1-4*l+3*l2)*it->A_w0 + (2*l-3*l2)*it->B_w0)*(dist  - 0.5*r_l - (r0))
					 -l*((1-2*l+l2)*it->A_w0 +l*(1-l)*it->B_w0)*(D_r0));
	  }

      }

    conf.current().perturbed_energy_derivatives.distrest_energy[topo.atom_energy_group()
								[it->v1.atom(0)]] += energy_derivativ;
    
  }

  return 0;
}

int interaction::Perturbed_Distance_Restraint_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{

  SPLIT_VIRIAL_BOUNDARY(_calculate_perturbed_distance_restraint_interactions,
			topo, conf, sim);
  
  return 0;
}
