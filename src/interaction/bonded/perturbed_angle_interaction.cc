/**
 * @file perturbed_angle_interaction.cc
 * template methods of Perturbed_Angle_Interaction
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
#include "perturbed_angle_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE bonded

/**
 * calculate angle forces and energies and lambda derivatives.
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_perturbed_angle_interactions
( topology::Topology & topo,
  configuration::Configuration & conf,
  simulation::Simulation & sim,
  interaction::Angle_Interaction const & m_interaction)
{
  // this is repeated code from Angle_Interaction !!!

  DEBUG(5, "perturbed angle interaction");
  DEBUG(7, "using the angle interaction: " << m_interaction.name);
  DEBUG(7, std::setprecision(5));
  
  // loop over the angles
  std::vector<topology::perturbed_three_body_term_struct>::const_iterator a_it =
    topo.perturbed_solute().angles().begin(),
    a_to = topo.perturbed_solute().angles().end();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec rij, rkj, fi, fj, fk;

  double energy, e_lambda;

  math::Periodicity<B> periodicity(conf.current().box);

  for( ; a_it != a_to; ++a_it){

    // atom i determines the energy group for the output. 
    // we use the same definition for the individual lambdas
    const double lambda = topo.individual_lambda(simulation::angle_lambda)
      [topo.atom_energy_group()[a_it->i]][topo.atom_energy_group()[a_it->i]];
    const double lambda_derivative = topo.individual_lambda_derivative
      (simulation::angle_lambda)
      [topo.atom_energy_group()[a_it->i]][topo.atom_energy_group()[a_it->i]];
    DEBUG(7, "angle " << a_it->i << "-" << a_it->j << "-" << a_it->k
	  << " A-type " << a_it->A_type
	  << " B-type " << a_it->B_type
	  << " lambda " << lambda);

    assert(pos.size() > (a_it->i) && pos.size() > (a_it->j) && 
	   pos.size() > (a_it->k));

    periodicity.nearest_image(pos(a_it->i), pos(a_it->j), rij);
    periodicity.nearest_image(pos(a_it->k), pos(a_it->j), rkj);

    double dij = sqrt(abs2(rij));
    double dkj = sqrt(abs2(rkj));
    
    assert(dij != 0.0);
    assert(dkj != 0.0);

    double ip = dot(rij, rkj);
    double cost = ip / (dij * dkj);

    assert(unsigned(a_it->A_type) < m_interaction.parameter().size());
    assert(unsigned(a_it->B_type) < m_interaction.parameter().size());
    
    double K    = (1 - lambda) *
      m_interaction.parameter()[a_it->A_type].K +
      lambda *
      m_interaction.parameter()[a_it->B_type].K;
    double cos0 =  (1 - lambda) *
      m_interaction.parameter()[a_it->A_type].cos0 +
      lambda *
      m_interaction.parameter()[a_it->B_type].cos0;

    const double K_diff = m_interaction.parameter()[a_it->B_type].K - 
      m_interaction.parameter()[a_it->A_type].K;
    const double cos_diff=m_interaction.parameter()[a_it->B_type].cos0- 
      m_interaction.parameter()[a_it->A_type].cos0;
    
    DEBUG(10, "K=" << K << " cos0=" << cos0 << " dij=" << dij << " dkj=" << dkj)
;

    double ki = -K * (cost - cos0) / dij;
    double kk = -K * (cost - cos0) / dkj;
    
    DEBUG(10, "cost=" << cost << " ki=" << ki << " kk=" << kk);

    fi = ki*(rkj/dkj - rij/dij * cost);
    fk = kk*(rij/dij - rkj/dkj * cost);
    fj = -1.0 * fi - fk;
    
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

    e_lambda = 0.5 * lambda_derivative * 
      ( -2.0 * K * cos_diff * (cost - cos0) +
	K_diff * (cost - cos0) * (cost - cos0) );

    DEBUG(9, "energy: " << energy);

    DEBUG(9, "K_diff: " << K_diff);

    DEBUG(9, "cos_diff: " << cos_diff);
    
    DEBUG(9, "e_lambda: " << e_lambda);
    
    assert(conf.current().energies.angle_energy.size() >
	   topo.atom_energy_group()[a_it->i]);
    
    conf.current().energies.
      angle_energy[topo.atom_energy_group()
		  [a_it->i]] += energy;
    
    assert(conf.current().perturbed_energy_derivatives.angle_energy.size() >
	   topo.atom_energy_group()[a_it->i]);
    
    conf.current().perturbed_energy_derivatives.
      angle_energy[topo.atom_energy_group()
		  [a_it->i]] += e_lambda;

  }

  return 0;
  
}

int interaction::Perturbed_Angle_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  m_timer.start();

  SPLIT_VIRIAL_BOUNDARY(_calculate_perturbed_angle_interactions,
			topo, conf, sim, m_interaction);

  m_timer.stop();

  return 0;
}
