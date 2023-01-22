/**
 * @file harm_angle_interaction.cc
 * template methods of Harm_Angle_Interaction.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"

#include "../../math/periodicity.h"

// interactions
#include "../../interaction/interaction_types.h"
#include "harm_angle_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE bonded

/**
 * calculate angle forces and energies.
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_harm_angle_interactions(topology::Topology & topo,
					      configuration::Configuration & conf,
					      simulation::Simulation & sim)
{
  std::vector<interaction::angle_type_struct> const & param = topo.angle_types_harm();
  // loop over the bonds
  std::vector<topology::three_body_term_struct>::const_iterator a_it =
    topo.solute().angles().begin(),
    a_to = topo.solute().angles().end();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec rij, rkj, fi, fj, fk;

  DEBUG(10, "pos   size: " << pos.size());
  DEBUG(10, "force size: " << force.size());
  
  double energy = 0.0;

  math::Periodicity<B> periodicity(conf.current().box);
  
  for( ; a_it != a_to; ++a_it){
    
    periodicity.nearest_image(pos(a_it->i), pos(a_it->j), rij);
    periodicity.nearest_image(pos(a_it->k), pos(a_it->j), rkj);

    DEBUG(9, "angle: " << a_it->i << " - " << a_it->j << " - " << a_it->k);

    DEBUG(10, "\trij = " << math::v2s(rij));
    DEBUG(10, "\trkj = " << math::v2s(rkj));

    double dij = sqrt(abs2(rij));
    double dkj = sqrt(abs2(rkj));

    DEBUG(10, "\tdij = " << dij);
    DEBUG(10, "\tdkj = " << dkj);
    
    assert(dij != 0.0);
    assert(dkj != 0.0);

    double ip = dot(rij, rkj);
    double cost = ip / (dij * dkj);

    DEBUG(10, "\tip = " << ip);

    double theta = acos(cost);
    double sint = sin(theta);
        
    assert(unsigned(a_it->type) < param.size());
 
    double K  = param[a_it->type].K;
    double theta0 = param[a_it->type].cos0;
    //here theta0 is actually cos0, as the parameter cosO is read differently depending on harmonicity 
    //see io::In_Topology::read_bondangle_types

    DEBUG(10, "\tK=" << K << " theta0=" << theta0 << " dij=" << dij << " dkj=" << dkj);

    double ki = 0.0, kk = 0.0;
    // double kj;

    if (sint < math::epsilon){
      // if sint is close to 0, which means theta is pi (as theta cannot be zero or -pi)
      if ((theta0 > math::Pi + math::epsilon) ||
	  (theta0 < math::Pi - math::epsilon)){
	io::messages.add("theta -> 180 but theta0 != 180",
			 "Harm_Angle_Interaction",
			 io::message::critical);
	// no force will be calculated
	continue;
      }
// to avoid numerical errors because (teta-teta0)/sinteta = 1 when teta=pi
      ki = -K / dij;
      // kj = -K;
      kk = -K / dkj;

      energy = K * (1 + cost);
    }
    else{
      ki = K * (theta - theta0) / sint / dij;
      // kj = K * (theta - theta0) / sint;
      kk = K * (theta - theta0) / sint / dkj;

      energy = 0.5 * K * (theta - theta0) * (theta - theta0);
    }
    
    fi = ki * (rkj / dkj - rij / dij * cost);
    fk = kk * (rij / dij - rkj / dkj * cost);
    fj = -fi - fk;
    // fj = kj * ((rij / (dij * dij) + rkj / (dkj * dkj)) * cost - (rij + rkj) / (dij * dkj));


    DEBUG(10, "\tenergy = " << energy);

    DEBUG(10, "\tfi=" << math::v2s(fi));
    DEBUG(10, "\tfj=" << math::v2s(fj));
    DEBUG(10, "\tfk=" << math::v2s(fk));
    
    force(a_it->i) += fi;
    force(a_it->j) += fj;
    force(a_it->k) += fk;

    // if (V == math::atomic_virial){
      for(int a=0; a<3; ++a)
	for(int bb=0; bb<3; ++bb)
	  conf.current().virial_tensor(a, bb) += 
	    rij(a) * fi(bb) +
	    rkj(a) * fk(bb);

      DEBUG(11, "\tatomic virial done");
      // }

    conf.current().energies.angle_energy[topo.atom_energy_group()[a_it->i]] += energy;
    // ORIOL_GAMD
    if(sim.param().gamd.gamd){
      unsigned int gamd_group = topo.gamd_accel_group(a_it->i);
      std::vector<unsigned int> key = {gamd_group, gamd_group};
      unsigned int igroup = topo.gamd_interaction_group(key);
      DEBUG(10, "\tGAMD interaction group is " << igroup);
      conf.special().gamd.total_force[igroup](a_it->i) += fi;
      conf.special().gamd.total_force[igroup](a_it->j) += fj;
      conf.special().gamd.total_force[igroup](a_it->k) += fk;
      conf.current().energies.gamd_potential_total[igroup] += energy;
      // virial
      for(int a=0; a<3; ++a){
        for(int bb=0; bb < 3; ++bb){
          conf.special().gamd.virial_tensor[igroup](a, bb) += rij(a) * fi(bb) + rkj(a) * fk(bb);
        }
      }

    } // end gamd

  }
  
  return 0;
}

int interaction::Harm_Angle_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  m_timer.start(sim);

  SPLIT_VIRIAL_BOUNDARY(_calculate_harm_angle_interactions,
			topo, conf, sim);

  m_timer.stop();

  return 0;
}
