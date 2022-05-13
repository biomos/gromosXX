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

  double energy = 0.0, e_lambda = 0.0;

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

    assert(unsigned(a_it->A_type) < topo.angle_types_cosharm().size());
    assert(unsigned(a_it->B_type) < topo.angle_types_cosharm().size());
    
    double K    = (1 - lambda) *
      topo.angle_types_cosharm()[a_it->A_type].K +
      lambda *
      topo.angle_types_cosharm()[a_it->B_type].K;
    double cos0 =  (1 - lambda) *
      topo.angle_types_cosharm()[a_it->A_type].cos0 +
      lambda *
      topo.angle_types_cosharm()[a_it->B_type].cos0;

    const double K_diff = topo.angle_types_cosharm()[a_it->B_type].K - 
      topo.angle_types_cosharm()[a_it->A_type].K;
    const double cos_diff=topo.angle_types_cosharm()[a_it->B_type].cos0- 
      topo.angle_types_cosharm()[a_it->A_type].cos0;
    
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

    // ANITA
    if (sim.param().precalclam.nr_lambdas &&
        ((sim.steps() % sim.param().write.free_energy) == 0)){
      double KA = topo.angle_types_cosharm()[a_it->A_type].K;
      double KB = topo.angle_types_cosharm()[a_it->B_type].K;
      double cos0A = topo.angle_types_cosharm()[a_it->A_type].cos0;
      double cos0B = topo.angle_types_cosharm()[a_it->B_type].cos0;

      double lambda_step = (sim.param().precalclam.max_lam -
                            sim.param().precalclam.min_lam) /
                            (sim.param().precalclam.nr_lambdas-1);

      //loop over nr_lambdas
      for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){

        // determine current lambda for this index
        double lam=(lam_index * lambda_step) + sim.param().precalclam.min_lam;

        double Klam = (1-lam)*KA + lam*KB;
        double cos0lam = (1-lam)*cos0A + lam*cos0B;
        double difflam = cost - cos0lam;
        double difflam2 = difflam * difflam;

        conf.current().energies.AB_angle[lam_index] += 0.5 * Klam * difflam2;
        conf.current().perturbed_energy_derivatives.AB_angle[lam_index] +=
                  0.5 * K_diff * difflam2 - Klam * cos_diff * difflam;
      }
    } //ANITA

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
