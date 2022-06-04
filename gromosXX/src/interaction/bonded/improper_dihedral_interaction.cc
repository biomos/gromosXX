/**
 * @file improper_dihedral_interaction.cc
 * template methods of Improper_dihedral_interaction.
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
#include "improper_dihedral_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE bonded

/**
 * calculate improper dihedral forces and energies.
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_improper_interactions(topology::Topology & topo,
					    configuration::Configuration & conf,
					    simulation::Simulation & sim)
{
  std::vector<interaction::improper_dihedral_type_struct> const & param = topo.impdihedral_types();
  // loop over the improper dihedrals
  std::vector<topology::four_body_term_struct>::const_iterator i_it =
    topo.solute().improper_dihedrals().begin(),
    i_to = topo.solute().improper_dihedrals().end();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec rij, rkj, rkl, rlj, rmj, rnk, fi, fj, fk, fl;
  double dkj2 = 0.0, dkj = 0.0, dmj2 = 0.0, dmj = 0.0, dnk2 = 0.0, dnk = 0.0, ip = 0.0, q = 0.0;
  double energy = 0.0;
  
  math::Periodicity<B> periodicity(conf.current().box);

  for( ; i_it != i_to; ++i_it){

    periodicity.nearest_image(pos(i_it->k), pos(i_it->j), rkj);
    periodicity.nearest_image(pos(i_it->i), pos(i_it->j), rij);
    periodicity.nearest_image(pos(i_it->k), pos(i_it->l), rkl);
    
    rmj = cross(rij, rkj);
    rnk = cross(rkj, rkl);
    
    dkj2 = abs2(rkj);
    dmj2 = abs2(rmj);
    dnk2 = abs2(rnk);
    dkj  = sqrt(dkj2);
    dmj  = sqrt(dmj2);
    dnk  = sqrt(dnk2);
    
    DEBUG(15,"dkj="<<dkj<<" dmj="<<dmj<<" dnk="<<dnk);
    
    assert(dmj != 0.0);
    assert(dnk != 0.0);

    ip = dot(rmj, rnk);
   
    double acs = ip / (dmj*dnk);
    if (acs > 1.0) {
      if (acs < 1.0 + math::epsilon) {
        acs = 1.0;
      } else {
        io::messages.add("improper dihedral",
                "acs > 1.0",
                io::message::critical);
      }
    }
    
    if (acs < -1.0) {
      if (acs > -1.0 - math::epsilon) {
        acs = -1.0;
      } else {
        io::messages.add("improper dihedral",
                "acs < -1.0",
                io::message::critical);
      }
    }
    
    q  = acos(acs);

    DEBUG(10, "zeta="<<q);
    
    ip = dot(rij, rnk);
    if(ip < 0) q *= -1.0;
    
    assert(unsigned(i_it->type) < param.size());
 
    const double K  = param[i_it->type].K;
    const double q0 = param[i_it->type].q0;

    double ki = -K * (q - q0) * dkj;
	double kl = -ki;
	if ( dmj2 < ( (1.0e-10 * dkj2))){
       ki = 0;
	   io::messages.add("One bond angle is close to 180 degrees!","improper_dihedral_interaction",io::message::warning);
    } else {
       ki = ki / dmj2;
    }
    if ( dnk2 < ( (1.0e-10 * dkj2))){
       kl = 0;
	   io::messages.add("One bond angle is close to 180 degrees!","improper_dihedral_interaction",io::message::warning);
    } else {
       kl = kl / dnk2;
    }

    const double kj1 = dot(rij, rkj) / dkj2 - 1.0;
    const double kj2 = dot(rkl, rkj) / dkj2;
    
    fi = ki * rmj;
    fl = kl * rnk;
    fj = kj1 * fi - kj2 * fl;
    fk = -1.0*(fi + fj + fl);
    
    force(i_it->i) += fi;
    force(i_it->j) += fj;
    force(i_it->k) += fk;
    force(i_it->l) += fl;
    
    // if (V == math::atomic_virial){
      periodicity.nearest_image(pos(i_it->l), pos(i_it->j), rlj);

      for(int a=0; a<3; ++a)
	for(int bb=0; bb<3; ++bb)
	  conf.current().virial_tensor(a, bb) += 
	    rij(a) * fi(bb) +
	    rkj(a) * fk(bb) +
	    rlj(a) * fl(bb);

      DEBUG(11, "\tatomic virial done");
      // }


    energy = 0.5 * K * (q-q0) * (q-q0);
    conf.current().energies.improper_energy[topo.atom_energy_group()[i_it->i]]
      += energy;
    
  }
  return 0;
  
}

int interaction::Improper_Dihedral_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  m_timer.start();
  
  SPLIT_VIRIAL_BOUNDARY(_calculate_improper_interactions,
			topo, conf, sim);

  m_timer.stop();

  return 0;
}
