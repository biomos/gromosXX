/**
 * @file perturbed_soft_improper_interaction.cc
 * template methods of Perturbed_Soft_Improper_Interaction
 */

#include <climits>
#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"

#include "../../math/periodicity.h"

// interactions
#include "../../interaction/interaction_types.h"
#include "perturbed_soft_improper_interaction.h"

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
static int _calculate_perturbed_soft_improper_interactions
( topology::Topology & topo,
  configuration::Configuration & conf,
  simulation::Simulation & sim)
{
  // this is repeated code from Improper_Dihedral_Interaction !!!
  std::vector<interaction::improper_dihedral_type_struct> &impropertypes = topo.impdihedral_types();

  DEBUG(5, "perturbed soft improper dihedral interaction");
  DEBUG(7, std::setprecision(5));
  
  // loop over the angles
  std::vector<topology::perturbed_four_body_term_struct>::const_iterator i_it =
    topo.perturbed_solute().softimpropers().begin(),
    i_to = topo.perturbed_solute().softimpropers().end();
    
  // and the softness parameters
  std::vector<double>::const_iterator alpha_it =
    topo.perturbed_solute().alpha_improper().begin();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec rij, rkj, rlj, rkl, rmj, rnk, fi, fj, fk, fl;
  double dkj2 = 0.0, dkj = 0.0, dmj2 = 0.0, dmj = 0.0, dnk2 = 0.0, dnk = 0.0, ip = 0.0, q = 0.0;
  double energy = 0.0, e_lambda = 0.0;

  math::Periodicity<B> periodicity(conf.current().box);

  for( ; i_it != i_to; ++i_it, ++alpha_it){

    // atom i determines the energy group for the output. 
    // we use the same definition for the individual lambdas
    const double lambda = topo.individual_lambda(simulation::improper_lambda)
      [topo.atom_energy_group()[i_it->i]][topo.atom_energy_group()[i_it->i]];
    const double lambda_derivative = topo.individual_lambda_derivative
      (simulation::improper_lambda)
      [topo.atom_energy_group()[i_it->i]][topo.atom_energy_group()[i_it->i]];

    DEBUG(7, "soft improper dihedral " << i_it->i << "-" << i_it->j << "-" 
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
    
    assert(unsigned(i_it->A_type) < impropertypes.size());

    const double K_A = impropertypes[i_it->A_type].K;
    const double K_B = impropertypes[i_it->B_type].K;
    double q0 =  (1 - lambda) *
      impropertypes[i_it->A_type].q0 +
      lambda *
      impropertypes[i_it->B_type].q0;
    double diff = q - q0;
    double diff2 = diff * diff;

    const double q_diff =
      impropertypes[i_it->B_type].q0- 
      impropertypes[i_it->A_type].q0;
    
    DEBUG(10, "K_A=" << K_A << ", K_B=" << K_B << ", q0=" << q0 );
    
    const double soft_A = 1 + *alpha_it * lambda * diff2;
    const double soft_B = 1 + *alpha_it * (1-lambda) * diff2;
    const double soft_A2 = soft_A*soft_A;
    const double soft_B2 = soft_B*soft_B;
    
    double fac = ((1-lambda)*K_A / soft_A2 + lambda*K_B / soft_B2) * diff * dkj;

    const double kj1 = dot(rij, rkj) / dkj2 - 1.0;
    const double kj2 = dot(rkl, rkj) / dkj2;
    
    fi = -fac / dmj2 * rmj;
    fl =  fac / dnk2 * rnk;
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

    energy = 0.5 * ((1-lambda)*K_A / soft_A + lambda*K_B / soft_B) * diff2;
    
    const double softterm1 = 1 + *alpha_it * diff2;
    const double softterm2 = -2 * *alpha_it * lambda * (1-lambda) * diff * q_diff;
    e_lambda = lambda_derivative
       * ( 0.5 * diff2 * ( K_A /  soft_A2  * ((-1) * softterm1 - softterm2) 
                        +  K_B /  soft_B2 * (softterm1 - softterm2))
       - diff * q_diff * ( K_A / soft_A * (1-lambda) + K_B / soft_B * lambda));
    
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

      double lambda_step = (sim.param().precalclam.max_lam -
                            sim.param().precalclam.min_lam) /
                            (sim.param().precalclam.nr_lambdas-1);

      //loop over nr_lambdas
      for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){

        // determine current lambda for this index
        double lam=(lam_index * lambda_step) + sim.param().precalclam.min_lam;

        const double q0lam = (1 - lam) * impropertypes[i_it->A_type].q0 
		                         + lam * impropertypes[i_it->B_type].q0;
        double difflam = q - q0lam;
        double difflam2 = difflam * difflam; 

        const double soft_Alam = 1 + *alpha_it * lam * difflam2;
        const double soft_Blam = 1 + *alpha_it * (1-lam) * difflam2;
        const double soft_Alam2 = soft_Alam*soft_Alam;
        const double soft_Blam2 = soft_Blam*soft_Blam;

        const double Ksoftlam = (1-lam)*K_A / soft_Alam + lam*K_B / soft_Blam;
        const double softterm1lam = 1 + *alpha_it * difflam2;
        const double softterm2lam = -2 * *alpha_it * lam * (1-lam) * difflam * q_diff;

        conf.current().energies.AB_improper[lam_index] += 0.5 * Ksoftlam * difflam2;
        conf.current().perturbed_energy_derivatives.AB_improper[lam_index] += 
                  0.5 * difflam2 * ( K_A /  soft_Alam2  * ((-1) * softterm1lam - softterm2lam) 
                        +  K_B /  soft_Blam2 * (softterm1lam - softterm2lam))
                  - Ksoftlam * difflam * q_diff;
      }
    } //ANITA
    
  }

  return 0;
  
}

int interaction::Perturbed_Soft_Improper_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  m_timer.start();
  
  SPLIT_VIRIAL_BOUNDARY(_calculate_perturbed_soft_improper_interactions,
			topo, conf, sim);

  m_timer.stop();

  return 0;
}



int interaction::Perturbed_Soft_Improper_Interaction
::init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os,
             bool quiet) 
    {
       if (!quiet)
           os << "Perturbed soft improper interaction\n";
      
      // add additional improper dih. types with K=0 and the target angle of 
      // the type we are perturbing to or from
      std::vector<topology::perturbed_four_body_term_struct>::iterator bt_it = topo.perturbed_solute().softimpropers().begin(),
                                 bt_to = topo.perturbed_solute().softimpropers().end();
      for( ; bt_it != bt_to; ++bt_it){
          if (bt_it->A_type==INT_MAX-1) {
            bt_it->A_type=topo.impdihedral_types().size();
            double qt = topo.impdihedral_types()[bt_it->B_type].q0;
            topo.impdihedral_types().push_back(interaction::improper_dihedral_type_struct(0, qt));
            DEBUG(10, "adding new improper dihedral type for soft improper perturbation: " 
                       << bt_it->A_type << " K=" << 0 << " qt=" << qt);                       
          } else if (bt_it->B_type==INT_MAX-1) {
            bt_it->B_type=topo.impdihedral_types().size();
            double qt = topo.impdihedral_types()[bt_it->A_type].q0;
            topo.impdihedral_types().push_back(interaction::improper_dihedral_type_struct(0, qt));   
            DEBUG(10, "adding new improper dihedral type for soft improper perturbation: " 
                       << bt_it->B_type << " K=" << 0 << " qt=" << qt);             
          }
      }
      return 0;
    }

