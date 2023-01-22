/**
 * @file angle_restraint_interaction.cc
 * template methods of Angle_Restraint_Interaction
 */

#include <sstream>
#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"

#include "../../math/periodicity.h"

// special interactions
#include "../../interaction/interaction_types.h"

#include "../../interaction/special/angle_restraint_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

/**
 * calculate angle restraint interactions
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_angle_restraint_interactions
(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim)
{
  DEBUG(5, "angle restraint interaction");
  // loop over the angle restraints
  std::vector<topology::angle_restraint_struct>::const_iterator 
    it = topo.angle_restraints().begin(),
    to = topo.angle_restraints().end();
    
   
  std::vector<double>::iterator ene_it = conf.special().angleres.energy.begin();
  std::vector<double>::iterator d_it = conf.special().angleres.d.begin();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec rij, rkj, fi, fj, fk;
  double energy = 0.0, f = 0.0;

  math::Periodicity<B> periodicity(conf.current().box);

  for(; it != to; ++it, ++ene_it, ++d_it){

    DEBUG(9, "angle " << it->i << "-" << it->j << "-" << it->k);

    periodicity.nearest_image(pos(it->k), pos(it->j), rkj);
    periodicity.nearest_image(pos(it->i), pos(it->j), rij);

    double dij = sqrt(abs2(rij));
    double dkj = sqrt(abs2(rkj));
    
    DEBUG(15,"dij="<<dij<<" dkj="<<dkj);

    assert(dij != 0.0);
    assert(dkj != 0.0);

    double ip = dot(rij, rkj);
    double cost = ip / (dij * dkj);
    
    if (cost > 1.0) {
      if (cost < 1.0 + math::epsilon) {
        cost = 1.0;
      } else {
        io::messages.add("angle",
                "cos(theta) > 1.0",
                io::message::critical);
      }
    }
    if (cost < -1.0) {
      if (cost > -1.0 - math::epsilon) {
        cost = -1.0;
      } else {
        io::messages.add("angle",
                "cos(theta) < -1.0",
                io::message::critical);
      }
    }
    
    double theta  = acos(cost);
    (*d_it) = theta;

    DEBUG(10, "raw theta="<< 180.0 * theta / math::Pi);
    
    double cost0 = it->cost;
    DEBUG(9, "theta=" << 180 * theta / math::Pi << " theta0=" << 180 * it->theta / math::Pi);

    double delta_theta = theta - it->theta;
    DEBUG(9, "delta_theta=" << 180 * delta_theta / math::Pi);
    double delta_cost = cost - cost0;

    double K = sim.param().angrest.K;
    if (sim.param().angrest.angrest == simulation::angle_restr_inst_weighted) {
      K *= it->w0;
    }
    DEBUG(10, "K=" << K);
    
    // HARMONIC
    DEBUG(10, "harmonic");
    energy = 0.5 * K * delta_cost * delta_cost;
    f = -K * delta_cost;
    

    DEBUG(10, "energy=" << energy << " force=" << f);
 
    (*ene_it) = energy;
    conf.current().energies.angrest_energy[topo.atom_energy_group()
					   [it->i]] += energy;

    fi = f / dij * (rkj/dkj - rij/dij * cost);
    fk = f / dkj * (rij/dij - rkj/dkj * cost);

    fj = -1.0 * fi - fk;

    DEBUG(10, "\tfi=" << math::v2s(fi));
    DEBUG(10, "\tfj=" << math::v2s(fj));
    DEBUG(10, "\tfk=" << math::v2s(fk));

    force(it->i) += fi;
    force(it->j) += fj;
    force(it->k) += fk;

    if (sim.param().angrest.virial) { 
      for (int a = 0; a < 3; ++a) {      
        for (int b = 0; b < 3; ++b) {
          conf.current().virial_tensor(a, b) +=
            rij(a) * fi(b) +
            rkj(a) * fk(b);
        }
      }
    }
    DEBUG(11, "\tatomic virial done");

  }
  
  return 0;
}

int interaction::Angle_Restraint_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{

  SPLIT_VIRIAL_BOUNDARY(_calculate_angle_restraint_interactions,
			topo, conf, sim);
  
  return 0;
}


  
 /**
 * initiate angle restraint interactions
 */
template<math::boundary_enum B>
static void _init_angres_data
(topology::Topology & topo,
 configuration::Configuration & conf)
{
   math::Periodicity<B> periodicity(conf.current().box);
   math::VArray &pos   = conf.current().pos;
 
   math::Vec rij, rkj;
    
  for(std::vector<topology::angle_restraint_struct>::const_iterator
        it = topo.angle_restraints().begin(),
        to = topo.angle_restraints().end(); it != to; ++it) {
        
    DEBUG(9, "init: angle " << it->i << "-" << it->j << "-" << it->k);

    
    periodicity.nearest_image(pos(it->i), pos(it->j), rij);
    periodicity.nearest_image(pos(it->k), pos(it->j), rkj);
      
    bool warn=false;
    for (int i=0; i<3;  i++) {
        if ((fabs(rij[i]) > conf.current().box(i)[i]*0.45 && abs(rij[i]) < conf.current().box(i)[i]*0.55)
         || (fabs(rkj[i]) > conf.current().box(i)[i]*0.45 && abs(rkj[i]) < conf.current().box(i)[i]*0.55)) {
          warn=true;
        }
    }
    if (warn) {
        std::ostringstream oss;
        oss << "one or more vectors of your angle restraint are\n"
          << "close to half a box length long in your initial structure, \n"
          << "the angle might be calculated from other periodic copies of the atoms than you intended!\n";
          io::messages.add(oss.str(),  "angle_restraint", io::message::warning);
    }

    double dij = sqrt(abs2(rij));
    double dkj = sqrt(abs2(rkj));
    
    DEBUG(15,"dij="<<dij<<" dkj="<<dkj);

    assert(dij != 0.0);
    assert(dkj != 0.0);

    double ip = dot(rij, rkj);
    double cost = ip / (dij * dkj);
    
    if (cost > 1.0) {
      if (cost < 1.0 + math::epsilon) {
        cost = 1.0;
      } else {
        io::messages.add("angle",
                "cos(theta) > 1.0",
                io::message::critical);
      }
    }
    if (cost < -1.0) {
      if (cost > -1.0 - math::epsilon) {
        cost = -1.0;
      } else {
        io::messages.add("angle",
                "cos(theta) < -1.0",
                io::message::critical);
      }
    }
    
    double theta  = acos(cost);
    
    DEBUG(10, "raw theta="<< 180.0 * theta / math::Pi);
             
    conf.special().angleres.d.push_back(theta);
    conf.special().angleres.energy.push_back(0.0);
  }
}

int interaction::Angle_Restraint_Interaction::init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os,
		     bool quiet) 
{
  
  SPLIT_BOUNDARY(_init_angres_data, topo, conf);
 
  if (!quiet) {
    os << "Angle restraint interaction";
    os << std::endl;
  }
  return 0;
}
