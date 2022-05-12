/**
 * @file perturbed_dihedral_new_interaction.cc
 * template methods of Perturbed_Dihedral_new_Interaction
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
#include "dihedral_new_interaction.h"
#include "perturbed_dihedral_new_interaction.h"

#include "../../util/template_split.h"
#include "../../util/error.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE bonded

#include "../../util/debug.h"

/**
 * calculate angle forces and energies and lambda derivatives.
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_perturbed_dihedral_new_interactions
(  topology::Topology & topo,
   configuration::Configuration & conf,
   simulation::Simulation & sim,
   interaction::Dihedral_new_Interaction const & m_interaction)
{
  // this is repeated code from Dihedral_new_Interaction !!!

  DEBUG(5, "perturbed dihedral interaction");
  DEBUG(7, "using the dihedral interaction: " 
	<< m_interaction.name);
  DEBUG(7, std::setprecision(5));
  
  // loop over the angles
  std::vector<topology::perturbed_four_body_term_struct>::const_iterator d_it =
    topo.perturbed_solute().dihedrals().begin(),
    d_to = topo.perturbed_solute().dihedrals().end();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec rij, rkj, rkl, rlj, rim, rln, rmj, rnk, fi, fj, fk, fl;
  math::Vec A_fi, A_fj, A_fk, A_fl, B_fi, B_fj, B_fk, B_fl;

  double dkj2 = 0.0, dkj = 0.0, dmj2 = 0.0, dnk2 = 0.0, dim = 0.0, dln = 0.0, ip = 0.0;
  
  double A_energy = 0.0, B_energy = 0.0, energy = 0.0, e_lambda = 0.0;

  math::Periodicity<B> periodicity(conf.current().box);

  for( ; d_it != d_to; ++d_it){

    // atom i determines the energy group for the output. 
    // we use the same definition for the individual lambdas
    const double lambda = topo.individual_lambda(simulation::dihedral_lambda)
      [topo.atom_energy_group()[d_it->i]][topo.atom_energy_group()[d_it->i]];
    const double lambda_derivative = topo.individual_lambda_derivative
      (simulation::dihedral_lambda)
      [topo.atom_energy_group()[d_it->i]][topo.atom_energy_group()[d_it->i]];
    
    DEBUG(7, "dihedral " << d_it->i << "-" << d_it->j << "-" 
	  << d_it->k << "-" << d_it->l
	  << " A-type " << d_it->A_type
	  << " B-type " << d_it->B_type
	  << " lambda " << lambda);
    
    assert(pos.size() > (d_it->i) && pos.size() > (d_it->j) && 
	   pos.size() > (d_it->k) && pos.size() > (d_it->l));

    periodicity.nearest_image(pos(d_it->k), pos(d_it->j), rkj);
    periodicity.nearest_image(pos(d_it->i), pos(d_it->j), rij);
    periodicity.nearest_image(pos(d_it->k), pos(d_it->l), rkl);

    rmj = cross(rij, rkj);
    rnk = cross(rkj, rkl);
    dmj2 = abs2(rmj);
    dnk2 = abs2(rnk);
    dkj2 = abs2(rkj);
    dkj = abs(rkj);
    

    double frim = dot(rij, rkj)/dkj2;
    double frln = dot(rkl, rkj)/dkj2;

    rim = rij - frim * rkj;
    rln = frln * rkj - rkl;
    dim = sqrt(abs2(rim));
    dln = sqrt(abs2(rln));
    
    ip = dot(rim, rln);

    double cosphi = ip / (dim*dln);
    if (cosphi > 1) cosphi = 1;
    if (cosphi < -1) cosphi = -1;
    double phi = acos(cosphi);
    double sign = dot(rij, rnk);
    if(sign < 0) phi*=-1.0;
    
    assert(unsigned(d_it->A_type) < topo.dihedral_types().size());
    
    // first state A 
   
    double     K = topo.dihedral_types()[d_it->A_type].K;
    double cosdelta = topo.dihedral_types()[d_it->A_type].cospd;
    double delta = topo.dihedral_types()[d_it->A_type].pd;
    double m = topo.dihedral_types()[d_it->A_type].m;
    
    DEBUG(10, "dihedral K=" << K << "cos delta=" << cosdelta);

    double ki = K * m * sin(m*phi - delta);
    double kl =  -ki;
    double kj1 = frim - 1.0;
    double kj2 = frln;
    
    A_fi = ki * dkj/dmj2 * rmj;
    A_fl = kl * dkj/dnk2 * rnk;
    A_fj = kj1 * A_fi - kj2 * A_fl;
    A_fk = -1.0 * (A_fi + A_fj + A_fl);
    
    A_energy = K * (1 + cos (m*phi - delta));
    
    // then state B 
    K = topo.dihedral_types()[d_it->B_type].K;
    delta = topo.dihedral_types()[d_it->B_type].pd;
    cosdelta = topo.dihedral_types()[d_it->B_type].cospd;
    m = topo.dihedral_types()[d_it->B_type].m;
    
    DEBUG(10, "dihedral K=" << K << "cos delta=" << cosdelta);

    ki = K * m * sin(m*phi - delta);
    kl =  -ki;
    kj1 = frim - 1.0;
    kj2 = frln;
    
    B_fi = ki * dkj/dmj2 * rmj;
    B_fl = kl * dkj/dnk2 * rnk;
    B_fj = kj1 * B_fi - kj2 * B_fl;
    B_fk = -1.0 * (B_fi + B_fj + B_fl);
    B_energy = K * (1 + cos (m*phi - delta));
    
    // now combine
    force(d_it->i) += (1.0 - lambda) * A_fi + lambda * B_fi;
    force(d_it->j) += (1.0 - lambda) * A_fj + lambda * B_fj;
    force(d_it->k) += (1.0 - lambda) * A_fk + lambda * B_fk;
    force(d_it->l) += (1.0 - lambda) * A_fl + lambda * B_fl;

    // if (V == math::atomic_virial){
      periodicity.nearest_image(pos(d_it->l), pos(d_it->j), rlj);

      for(int a=0; a<3; ++a)
	for(int bb=0; bb<3; ++bb)
	  conf.current().virial_tensor(a, bb) += 
	    rij(a) * ((1.0-lambda) * A_fi(bb) + lambda * B_fi(bb)) +
	    rkj(a) * ((1.0-lambda) * A_fk(bb) + lambda * B_fk(bb)) +
	    rlj(a) * ((1.0-lambda) * A_fl(bb) + lambda * B_fl(bb));

      DEBUG(11, "\tatomic virial done");
      // }

    energy = (1.0 - lambda) * A_energy + lambda * B_energy;
    e_lambda = lambda_derivative * (B_energy - A_energy);
    DEBUG(10, "energy " << energy << " e_lambda " << e_lambda);
    DEBUG(10, "force i " << force(d_it->i)[0] << " " << force(d_it->i)[1] << " " << force(d_it->i)[2]);
    DEBUG(10, "force j " << force(d_it->j)[0] << " " << force(d_it->j)[1] << " " << force(d_it->j)[2]);
    DEBUG(10, "force k " << force(d_it->k)[0] << " " << force(d_it->k)[1] << " " << force(d_it->k)[2]);
    DEBUG(10, "force l " << force(d_it->l)[0] << " " << force(d_it->l)[1] << " " << force(d_it->l)[2]);


    assert(conf.current().energies.dihedral_energy.size() >
	   topo.atom_energy_group()[d_it->i]);
    conf.current().energies.dihedral_energy
      [topo.atom_energy_group()[d_it->i]] += energy;
    
    assert(conf.current().perturbed_energy_derivatives.dihedral_energy.size() >
	   topo.atom_energy_group()[d_it->i]);
    conf.current().perturbed_energy_derivatives.dihedral_energy
      [topo.atom_energy_group()[d_it->i]] += e_lambda;   

    // ANITA
    if (sim.param().precalclam.nr_lambdas &&
        ((sim.steps() % sim.param().write.free_energy) == 0)){

      conf.current().energies.A_dihedral += A_energy;
      conf.current().energies.B_dihedral += B_energy;

      conf.current().perturbed_energy_derivatives.A_dihedral += B_energy - A_energy;


/*      double lambda_step = (sim.param().precalclam.max_lam -
                            sim.param().precalclam.min_lam) /
                            (sim.param().precalclam.nr_lambdas-1);

      //loop over nr_lambdas
      for (int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){

        // determine current lambda for this index
        double lam=(lam_index * lambda_step) + sim.param().precalclam.min_lam;

        conf.current().energies.AB_dihedral[lam_index] += (1-lam)*A_energy + lam*B_energy;
        conf.current().perturbed_energy_derivatives.AB_dihedral[lam_index] +=
                  B_energy - A_energy;
      } */
    } //ANITA

  }

  return 0;

}

int interaction::Perturbed_Dihedral_new_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation & sim)
{
  m_timer.start();
  
  SPLIT_VIRIAL_BOUNDARY(_calculate_perturbed_dihedral_new_interactions,
			topo, conf, sim, m_interaction);

  m_timer.stop();

  return 0;
}
