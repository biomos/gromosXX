/**
 * @file distance_restraint_interaction.cc
 * template methods of Distance_Restraint_Interaction
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

#include "../../interaction/special/distance_restraint_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

/**
 * calculate distance restraint interactions
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_distance_restraint_interactions
(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim, double exponential_term,
std::map<int,math::Vec> &rah_map, int &error)
{
  // loop over the distance restraints
  std::vector<topology::distance_restraint_struct>::const_iterator 
    it = topo.distance_restraints().begin(),
    to = topo.distance_restraints().end();
  
 std::vector<double>::iterator ave_it = conf.special().distanceres.av.begin();
 std::vector<double>::iterator ene_it = conf.special().distanceres.energy.begin();
 std::vector<double>::iterator d_it = conf.special().distanceres.d.begin();

  // math::VArray &pos   = conf.current().pos;
  // math::VArray &force = conf.current().force;
  math::Vec v, f;
  double energy = 0.0;

  math::Periodicity<B> periodicity(conf.current().box);

  // for(int i=0; it != to; ++it, ++i){
  for(; it != to; ++it, ++ave_it, ++ene_it, ++d_it) {
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
    DEBUG(9, "DISTANCERES rah: " << it->rah);

    // determine the dimensionality and a local copy of rah
    int rah=it->rah;
    bool found=false;
    for(std::map<int,math::Vec>::iterator i = rah_map.begin(), to = rah_map.end(); i != to; i++){
      if(it->rah>=(i->first) - 1 && it->rah <= (i->first) + 1){
        for(int j=0; j<3; j++){
          v[j] = v[j]*(i->second)[j];
        }
        rah = it->rah - (i->first);
        found = true;
      }
    }
    if(!found){
       std::ostringstream msg;
       msg << "Do not know how to handle " << it->rah << " for RAH in distance restraint" << std::endl;
       io::messages.add(msg.str(),
                         "Distance_restraints",
                         io::message::error);
       error = 1;
       return error;
    }
    DEBUG(9, "DISTANCERES updated rah: " << rah);
    DEBUG(9, "DISTANCERES updated v : " << math::v2s(v));

    double dist = math::abs(v);
    // dist is possibly overwritten with the time-averaged distance
    // store the instanteneous value in vabs
    double vabs = dist;
    
    (*d_it) = dist;
    DEBUG(9, "DISTANCERES dist : " << dist << " r0 : " << it->r0);
    double force_scale = 1.0;
    if (sim.param().distanceres.distanceres < 0) {

      (*ave_it) = (1.0 - exponential_term) * pow(dist, -3.0) + 
                   exponential_term * (*ave_it);

      dist = pow(*ave_it, -1.0 / 3.0);
      DEBUG(9, "DISTANCERES average dist : " << dist);

      // the force_scaling causes large fluctuations and may be omitted
      if(sim.param().distanceres.forcescale == 0)
	force_scale = 1.0;
      else if(sim.param().distanceres.forcescale == 1)
	force_scale = (1.0 - exponential_term);
      else if(sim.param().distanceres.forcescale == 2)
	force_scale = (1.0 - exponential_term) * dist * dist * dist * dist / 
	  (vabs * vabs * vabs * vabs);
    } else {
      (*ave_it) = 0.0;
    }

    if(rah*dist < rah * it->r0){
      DEBUG(9, "DISTANCERES  : restraint fulfilled");
      energy = 0;
      f=0;
    }
    else if(fabs(it->r0 - dist) < sim.param().distanceres.r_linear){
      DEBUG(9, "DISTANCERES  : harmonic");
      energy = 0.5 * sim.param().distanceres.K * (dist - it->r0) * (dist - it->r0);
      f = - sim.param().distanceres.K * (dist - it->r0) * v / vabs;
    }
    else{
      DEBUG(9, "DISTANCERES  : linear");
      if(dist < it->r0){
	energy = -sim.param().distanceres.K * sim.param().distanceres.r_linear *
	  (dist - it->r0 + 0.5 * sim.param().distanceres.r_linear);
	f =  sim.param().distanceres.r_linear * sim.param().distanceres.K  * v / vabs;
      }
      else {
	energy =  sim.param().distanceres.K * sim.param().distanceres.r_linear *
	  (dist - it->r0 - 0.5 * sim.param().distanceres.r_linear); 
	f = -sim.param().distanceres.r_linear * sim.param().distanceres.K  * v / vabs;
      }
    }
    
    if(abs(sim.param().distanceres.distanceres) == 1)
      ;
    else if(abs(sim.param().distanceres.distanceres) == 2){
      energy *= it->w0;
      f *= it->w0;
    }
    else{
      energy = 0;
      f = 0;
    }
    
    // scale the force according to 8.17
    f *= force_scale;
    DEBUG(9, "Distanceres force : " << math::v2s(f));

    it->v1.force(conf, topo,  f);
    it->v2.force(conf, topo, -f);  

    (*ene_it) = energy;

    DEBUG(9, "Distanceres energy : " << energy);

    conf.current().energies.distanceres_energy[topo.atom_energy_group()
					    [it->v1.atom(0)]] += energy;

    if (sim.param().distanceres.virial) { 
      for (int a = 0; a < 3; ++a) {
	for (int b = 0; b < 3; ++b) { 
	  conf.current().virial_tensor(a, b) += v(a) * f(b); 
	}
      } 
    }
       
  }
  
  return 0;
}

int interaction::Distance_Restraint_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  int error =0;
  SPLIT_VIRIAL_BOUNDARY(_calculate_distance_restraint_interactions,
			topo, conf, sim, exponential_term, rah_map, error);
  
  return error;
}

/**
 * initiate distance restraint interactions
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
    periodicity.nearest_image(it->v1.pos(conf,topo), it->v2.pos(conf,topo), v);
    conf.special().distanceres.d.push_back(math::abs(v));
    conf.special().distanceres.energy.push_back(0.0);
    conf.special().distanceres.av.push_back(pow(math::abs(v), -3.0));
  }
}

/**
 * initiate distance restraint interactions
 */
template<math::boundary_enum B>
static void _init_disres_data
(topology::Topology & topo,
 configuration::Configuration & conf)
{
  math::Periodicity<B> periodicity(conf.current().box);
  math::Vec v;

  for(std::vector<topology::distance_restraint_struct>::const_iterator
        it = topo.distance_restraints().begin(),
        to = topo.distance_restraints().end(); it != to; ++it) {
    periodicity.nearest_image(it->v1.pos(conf,topo), it->v2.pos(conf,topo), v);
    conf.special().distanceres.d.push_back(math::abs(v));
    conf.special().distanceres.energy.push_back(0.0);
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
    } else {
      SPLIT_BOUNDARY(_init_disres_data, topo, conf);
    }
  } else {
    SPLIT_BOUNDARY(_init_averages, topo, conf);
  }
 
  // set the dimensionality map only once 
  rah_map[0] = math::Vec(1,1,1);
  rah_map[10] = math::Vec(1,1,0);
  rah_map[20] = math::Vec(1,0,1);
  rah_map[30] = math::Vec(0,1,1);
  rah_map[40] = math::Vec(1,0,0);
  rah_map[50] = math::Vec(0,1,0);
  rah_map[60] = math::Vec(0,0,1);
 
  if (!quiet) {
    os << "Distance restraint interaction";
    if (sim.param().distanceres.distanceres < 0)
      os << " with time-averaging";
    os << std::endl;
  }
  return 0;
}
