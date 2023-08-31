/**
 * @file angle_interaction.cc
 * template methods of Angle_Interaction.
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
#include "angle_interaction.h"

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
static int _calculate_angle_interactions(topology::Topology & topo,
					 configuration::Configuration & conf,
					 simulation::Simulation & sim)
{
  std::vector<interaction::angle_type_struct> const & param = topo.angle_types_cosharm();
  // loop over the bonds
  std::vector<topology::three_body_term_struct>::const_iterator
    a_it = topo.solute().angles().begin(),
    a_to = topo.solute().angles().end();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec rij, rkj, fi, fj, fk;

  DEBUG(10, "pos   size: " << pos.size());
  DEBUG(10, "force size: " << force.size());
  
  double energy = 0.0;

  DEBUG(7, "bond angles : " << topo.solute().angles().size());

  math::Periodicity<B> periodicity(conf.current().box);
  
  for( ; a_it != a_to; ++a_it){
    
    periodicity.nearest_image(pos(a_it->i), pos(a_it->j), rij);
    periodicity.nearest_image(pos(a_it->k), pos(a_it->j), rkj);

    DEBUG(9, "g96 angle: " << a_it->i << " - " << a_it->j << " - " << a_it->k);

    DEBUG(10, "\tpos(i) = " << math::v2s(pos(a_it->i)));
    DEBUG(10, "\tpos(j) = " << math::v2s(pos(a_it->j)));
    DEBUG(10, "\tpos(k) = " << math::v2s(pos(a_it->k)));

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
        
    assert(unsigned(a_it->type) < param.size());
 
    double K    = param[a_it->type].K;
    double cos0 = param[a_it->type].cos0;

    DEBUG(10, "\tK=" << K << " cos0=" << cos0 << " dij=" << dij << " dkj=" << dkj);

    const double df = -K * (cost - cos0);

    DEBUG(10, "\tcost=" << cost << " df=" << df);

    fi = df / dij * (rkj/dkj - rij/dij * cost);
    fk = df / dkj * (rij/dij - rkj/dkj * cost);

    fj = -1.0 * fi - fk;

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


    energy = 0.5 * K * (cost - cos0) * (cost - cos0);
    conf.current().energies.angle_energy[topo.atom_energy_group()[a_it->i]]
      += energy;

    DEBUG(10, "\tenergy = " << energy);
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

int interaction::Angle_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  int e=0;
  m_timer.start(sim);

  SPLIT_VIRIAL_BOUNDARY(_calculate_angle_interactions,
			topo, conf, sim);

  m_timer.stop();

  return e;
}
