/**
 * @file perturbed_dihedral_interaction.cc
 * template methods of Perturbed_Dihedral_Interaction
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
#include "dihedral_interaction.h"
#include "perturbed_dihedral_interaction.h"

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
static int _calculate_perturbed_dihedral_interactions
(  topology::Topology & topo,
   configuration::Configuration & conf,
   simulation::Simulation & sim,
   interaction::Dihedral_Interaction const & m_interaction)
{
  // this is repeated code from Dihedral_Interaction !!!

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

  double dkj2 = 0.0, dim = 0.0, dln = 0.0, ip = 0.0;
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
    
    dkj2 = abs2(rkj);

    double frim = dot(rij, rkj)/dkj2;
    double frln = dot(rkl, rkj)/dkj2;

    rim = rij - frim * rkj;
    rln = frln * rkj - rkl;
    dim = sqrt(abs2(rim));
    dln = sqrt(abs2(rln));
    
    ip = dot(rim, rln);

    double cosphi = ip / (dim*dln);
    
    double cosphi2 = cosphi  * cosphi;
    double cosphi3 = cosphi2 * cosphi;
    double cosphi4 = cosphi3 * cosphi;

    assert(unsigned(d_it->A_type) < topo.dihedral_types().size());
    
    double dcosmphi = 0;
    double cosmphi = 0;

    // first state A 
    switch(topo.dihedral_types()[d_it->A_type].m){
      case 0:
        cosmphi = 0.0;
        dcosmphi = 0.0;
        break;
      case 1:
        cosmphi = cosphi;
        dcosmphi = 1;
        break;
      case 2:
        cosmphi =  2*cosphi2 -1;
        dcosmphi = 4*cosphi;
        break;
      case 3:
        cosmphi  = 4*cosphi3 - 3*cosphi;
        dcosmphi = 12*cosphi2 - 3;
        break;
      case 4:
        cosmphi  = 8*cosphi4 - 8*cosphi2 + 1;
        dcosmphi = 32*cosphi3-16*cosphi;
        break;
      case 5:
        cosmphi  = 16*cosphi4*cosphi - 20*cosphi3 + 5*cosphi;
        dcosmphi = 80*cosphi4-60*cosphi2+5;
        break;
      case 6:
        cosmphi  = 32*cosphi4*cosphi2 - 48*cosphi4 + 18*cosphi2 -1;
        dcosmphi = 192*cosphi4*cosphi-192*cosphi3+36*cosphi;
        break;
      
    }
    double     K = topo.dihedral_types()[d_it->A_type].K;
    double     cosdelta = topo.dihedral_types()[d_it->A_type].cospd;
    
    DEBUG(10, "dihedral K=" << K << "cos delta=" << cosdelta << " dcos=" << dcosmphi);

    double ki = -K * cosdelta * dcosmphi / dim;
    double kl = -K * cosdelta * dcosmphi / dln;
    double kj1 = frim - 1.0;
    double kj2 = frln;
    
    A_fi = ki * (rln / dln - rim / dim * cosphi);
    A_fl = kl * (rim / dim - rln / dln * cosphi);
    A_fj = kj1 * A_fi - kj2 * A_fl;
    A_fk = -1.0 * (A_fi + A_fj + A_fl);
    

    A_energy = K * (1 + cosdelta * cosmphi);

    // then state B 
    switch(topo.dihedral_types()[d_it->B_type].m){
      case 0:
        cosmphi = 0.0;
        dcosmphi = 0.0;
        break;
      case 1:
        cosmphi = cosphi;
        dcosmphi = 1;
        break;
      case 2:
        cosmphi =  2*cosphi2 -1;
        dcosmphi = 4*cosphi;
        break;
      case 3:
        cosmphi  = 4*cosphi3 - 3*cosphi;
        dcosmphi = 12*cosphi2 - 3;
        break;
      case 4:
        cosmphi  = 8*cosphi4 - 8*cosphi2 + 1;
        dcosmphi = 32*cosphi3-16*cosphi;
        break;
      case 5:
        cosmphi  = 16*cosphi4*cosphi - 20*cosphi3 + 5*cosphi;
        dcosmphi = 80*cosphi4-60*cosphi2+5;
        break;
      case 6:
        cosmphi  = 32*cosphi4*cosphi2 - 48*cosphi4 + 18*cosphi2 -1;
        dcosmphi = 192*cosphi4*cosphi-192*cosphi3+36*cosphi;
        break;
       }
        K = topo.dihedral_types()[d_it->B_type].K;
    cosdelta = topo.dihedral_types()[d_it->B_type].cospd;

    DEBUG(10, "dihedral K=" << K << "cos delta=" << cosdelta << " dcos=" << dcosmphi);

    ki = -K * cosdelta * dcosmphi / dim;
    kl = -K * cosdelta * dcosmphi / dln;
    kj1 = frim - 1.0;
    kj2 = frln;
    
    B_fi = ki * (rln / dln - rim / dim * cosphi);
    B_fl = kl * (rim / dim - rln / dln * cosphi);
    B_fj = kj1 * B_fi - kj2 * B_fl;
    B_fk = -1.0 * (B_fi + B_fj + B_fl);
    
    B_energy = K * (1 + cosdelta * cosmphi);

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

int interaction::Perturbed_Dihedral_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation & sim)
{
  m_timer.start();
  
  SPLIT_VIRIAL_BOUNDARY(_calculate_perturbed_dihedral_interactions,
			topo, conf, sim, m_interaction);

  m_timer.stop();

  return 0;
}

int interaction::Perturbed_Dihedral_Interaction::init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os,
                     bool quiet) {
       
    std::vector<interaction::dihedral_type_struct> const & param = topo.dihedral_types();
    std::vector<topology::four_body_term_struct>::iterator d_it =
    topo.solute().dihedrals().begin(),
    d_to = topo.solute().dihedrals().end();
  
    for(int n =0; d_it != d_to; ++d_it, ++n){
        double cosdelta = param[d_it->type].cospd;
        int m = param[d_it->type].m;
  
       if (((cosdelta > -1 - math::epsilon)  &&  (cosdelta < -1 + math::epsilon)) || ((cosdelta > 1 - math::epsilon)  &&  (cosdelta < 1 + math::epsilon))) {
        //          
       }
       else {
            io::messages.add("perturbed dihedral function (NTBDN=0) not implemented for phase shifts not equal to 0 or 180 degrees. Please use NTBDN=1 in COVALENTFORM block", "perturbed_dihedral interaction", io::message::error);
            return 1;
            } 
        if (m>6 || m<0) {
            io::messages.add("perturbed dihedral function not implemented for m>6 or m<0. Please use NTBDN=1 in COVALENTFORM block", "perturbed_dihedral_interaction", io::message::error);
            return 1;
        }
        }
    return 0;

 };
