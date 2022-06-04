/**
 * @file perturbed_improper_dihedral_interaction.cc
 * template methods of Perturbed_Improper_Dihedral_Interaction
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
#include "perturbed_improper_dihedral_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE bonded

/**
 * calculate improper dihedral forces and energies and lambda derivatives.
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_perturbed_improper_interactions
( topology::Topology & topo,
  configuration::Configuration & conf,
  simulation::Simulation & sim,
  interaction::Improper_Dihedral_Interaction const & m_interaction)
{
  // this is repeated code from Improper_Dihedral_Interaction !!!

  DEBUG(5, "perturbed improper dihedral interaction");
  DEBUG(7, "using the improper dihedral interaction: " 
	<< m_interaction.name);
  DEBUG(7, std::setprecision(5));
  
  // loop over the angles
  std::vector<topology::perturbed_four_body_term_struct>::const_iterator i_it =
    topo.perturbed_solute().improper_dihedrals().begin(),
    i_to = topo.perturbed_solute().improper_dihedrals().end();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec rij, rkj, rlj, rkl, rmj, rnk, fi, fj, fk, fl;
  double dkj2 = 0.0, dkj = 0.0, dmj2 = 0.0, dmj = 0.0, dnk2 = 0.0, dnk = 0.0, ip = 0.0, q = 0.0;
  double energy = 0.0, e_lambda = 0.0;

  math::Periodicity<B> periodicity(conf.current().box);

  for( ; i_it != i_to; ++i_it){

    // atom i determines the energy group for the output. 
    // we use the same definition for the individual lambdas
    const double lambda = topo.individual_lambda(simulation::improper_lambda)
      [topo.atom_energy_group()[i_it->i]][topo.atom_energy_group()[i_it->i]];
    const double lambda_derivative = topo.individual_lambda_derivative
      (simulation::improper_lambda)
      [topo.atom_energy_group()[i_it->i]][topo.atom_energy_group()[i_it->i]];

    DEBUG(7, "improper dihedral " << i_it->i << "-" << i_it->j << "-" 
	  << i_it->k << "-" << i_it->l
	  << " A-type " << i_it->A_type
	  << " B-type " << i_it->B_type
	  << " lambda " << lambda);

    assert(pos.size() > (i_it->i) && pos.size() > (i_it->j) && 
	   pos.size() > (i_it->k) && pos.size() > (i_it->l));

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
   
    q  = acos(ip / (dmj*dnk));

    DEBUG(10, "zeta="<<q);
    
    ip = dot(rij, rnk);
    if(ip < 0) q *= -1.0;
    
    assert(unsigned(i_it->A_type) < topo.impdihedral_types().size());

    double K    = (1 - lambda) *
      topo.impdihedral_types()[i_it->A_type].K +
      lambda *
      topo.impdihedral_types()[i_it->B_type].K;
    double q0 =  (1 - lambda) *
      topo.impdihedral_types()[i_it->A_type].q0 +
      lambda *
      topo.impdihedral_types()[i_it->B_type].q0;

    const double K_diff = 
      topo.impdihedral_types()[i_it->B_type].K - 
      topo.impdihedral_types()[i_it->A_type].K;
    const double q_diff =
      topo.impdihedral_types()[i_it->B_type].q0- 
      topo.impdihedral_types()[i_it->A_type].q0;
    
    DEBUG(10, "K=" << K << " q0=" << q0 );

    const double ki = -K * (q - q0) * dkj / dmj2;
    const double kl = K * (q - q0) * dkj / dnk2;
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
    
    e_lambda = 0.5 * lambda_derivative* ( -2.0 * K * q_diff * (q-q0) +
					  K_diff * (q-q0) * (q-q0));
    
    assert(conf.current().energies.improper_energy.size() >
	   topo.atom_energy_group()[i_it->i]);
    conf.current().energies.
      improper_energy[topo.atom_energy_group()
		      [i_it->i]] += energy;

    assert(conf.current().perturbed_energy_derivatives.improper_energy.size() >
	   topo.atom_energy_group()[i_it->i]);
    
    conf.current().perturbed_energy_derivatives.
      improper_energy[topo.atom_energy_group()
		      [i_it->i]] += e_lambda;

    // ANITA
    if (sim.param().precalclam.nr_lambdas &&
        ((sim.steps() % sim.param().write.free_energy) == 0)){
      double KA = topo.impdihedral_types()[i_it->A_type].K;
      double KB = topo.impdihedral_types()[i_it->B_type].K;
      double q0A = topo.impdihedral_types()[i_it->A_type].q0;
      double q0B = topo.impdihedral_types()[i_it->B_type].q0;

      double lambda_step = (sim.param().precalclam.max_lam -
                            sim.param().precalclam.min_lam) /
                            (sim.param().precalclam.nr_lambdas-1);

      //loop over nr_lambdas
      for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){

        // determine current lambda for this index
        double lam=(lam_index * lambda_step) + sim.param().precalclam.min_lam;

        double Klam = (1-lam)*KA + lam*KB;
        double q0lam = (1-lam)*q0A + lam*q0B;
        double difflam = q - q0lam;
        double difflam2 = difflam * difflam;

        conf.current().energies.AB_improper[lam_index] += 0.5 * Klam * difflam2;
        conf.current().perturbed_energy_derivatives.AB_improper[lam_index] +=
                  0.5 * K_diff * difflam2 - Klam * q_diff * difflam;
      }
    } //ANITA
    
  }

  return 0;
  
}

int interaction::Perturbed_Improper_Dihedral_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  m_timer.start();
  
  SPLIT_VIRIAL_BOUNDARY(_calculate_perturbed_improper_interactions,
			topo, conf, sim, m_interaction);

  m_timer.stop();

  return 0;
}
