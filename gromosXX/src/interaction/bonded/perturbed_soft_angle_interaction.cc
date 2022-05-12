/**
 * @file perturbed_soft_angle_interaction.cc
 * template methods of Perturbed_Soft_Angle_Interaction
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
#include "perturbed_soft_angle_interaction.h"

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
static int _calculate_perturbed_soft_angle_interactions
( topology::Topology & topo,
  configuration::Configuration & conf,
  simulation::Simulation & sim)
{
  // this is repeated code from Angle_Interaction !!!
  std::vector<interaction::angle_type_struct> &angletypes = topo.angle_types_cosharm();

  DEBUG(5, "perturbed soft angle interaction");
  DEBUG(7, std::setprecision(5));
  
  // loop over the angles
  std::vector<topology::perturbed_three_body_term_struct>::const_iterator a_it =
    topo.perturbed_solute().softangles().begin(),
    a_to = topo.perturbed_solute().softangles().end();
    
  // and the softness parameters
  std::vector<double>::const_iterator alpha_it =
    topo.perturbed_solute().alpha_angle().begin();


  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec rij, rkj, fi, fj, fk;

  double energy = 0.0, e_lambda = 0.0;

  math::Periodicity<B> periodicity(conf.current().box);

  for( ; a_it != a_to; ++a_it,++alpha_it){

    // atom i determines the energy group for the output. 
    // we use the same definition for the individual lambdas
    const double lambda = topo.individual_lambda(simulation::angle_lambda)
      [topo.atom_energy_group()[a_it->i]][topo.atom_energy_group()[a_it->i]];
    const double lambda_derivative = topo.individual_lambda_derivative
      (simulation::angle_lambda)
      [topo.atom_energy_group()[a_it->i]][topo.atom_energy_group()[a_it->i]];
    DEBUG(7, "soft angle " << a_it->i << "-" << a_it->j << "-" << a_it->k
	  << " A-type " << a_it->A_type
	  << " B-type " << a_it->B_type
	  << " alpha " << *alpha_it
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

    assert(unsigned(a_it->A_type) < angletypes.size());
    assert(unsigned(a_it->B_type) < angletypes.size());
    
    const double K_A = angletypes[a_it->A_type].K;
    const double K_B = angletypes[a_it->B_type].K;
    
    double cos0 =  (1 - lambda) *
      angletypes[a_it->A_type].cos0 +
      lambda *
      angletypes[a_it->B_type].cos0;
    double diff = cost - cos0;
    double diff2 = diff * diff;

    const double K_diff = K_B-K_A;
    const double cos_diff=angletypes[a_it->B_type].cos0- 
      angletypes[a_it->A_type].cos0;
    
    DEBUG(10, "K_A=" << K_A << " K_B=" << K_B << " cos0=" << cos0 << " dij=" << dij << " dkj=" << dkj);
    
    const double soft_A = 1 + *alpha_it * lambda * diff2;
    const double soft_B = 1 + *alpha_it * (1-lambda) * diff2;
    const double soft_A2 = soft_A*soft_A;
    const double soft_B2 = soft_B*soft_B;
    
    
    DEBUG(10, "cost=" << cost << " soft_A=" << soft_A << " soft_B=" << soft_B);

    double fac = ((1-lambda)*K_A / soft_A2 + lambda*K_B / soft_B2) * diff;
    
    fi = -fac * (rkj/dkj - rij/dij * cost) / dij;
    fk = -fac * (rij/dij - rkj/dkj * cost) / dkj;
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
    const double Ksoft = (1-lambda)*K_A / soft_A + lambda*K_B / soft_B;
    energy = 0.5 * Ksoft * diff2;

    const double softterm1 = 1 + *alpha_it * diff2;
    const double softterm2 = -2 * *alpha_it * lambda * (1-lambda) * diff * cos_diff;
    
    e_lambda = lambda_derivative 
       * ( 0.5 * diff2 * ( K_A /  soft_A2  * ((-1) * softterm1 - softterm2) 
                        +  K_B /  soft_B2 * (softterm1 - softterm2))
       - diff * cos_diff * Ksoft);

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

      double lambda_step = (sim.param().precalclam.max_lam -
                            sim.param().precalclam.min_lam) /
                            (sim.param().precalclam.nr_lambdas-1);

      //loop over nr_lambdas
      for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){

        // determine current lambda for this index
        double lam=(lam_index * lambda_step) + sim.param().precalclam.min_lam;

        double cos0lam = (1-lam)*angletypes[a_it->A_type].cos0
                           + lam*angletypes[a_it->B_type].cos0;
        double difflam = cost - cos0lam;
        double difflam2 = difflam * difflam;

        const double soft_Alam = 1 + *alpha_it * lam * difflam2;
        const double soft_Blam = 1 + *alpha_it * (1-lam) * difflam2;
        const double soft_Alam2 = soft_Alam*soft_Alam;
        const double soft_Blam2 = soft_Blam*soft_Blam;

        const double Ksoftlam = (1-lam)*K_A / soft_Alam + lam*K_B / soft_Blam;
        const double softterm1lam = 1 + *alpha_it * difflam2;
        const double softterm2lam = -2 * *alpha_it * lam * (1-lam) * difflam * cos_diff;

        conf.current().energies.AB_angle[lam_index] += 0.5 * Ksoftlam * difflam2;
        conf.current().perturbed_energy_derivatives.AB_angle[lam_index] +=
                  0.5 * difflam2 * ( K_A /  soft_Alam2  * ((-1) * softterm1lam - softterm2lam) 
                        +  K_B /  soft_Blam2 * (softterm1lam - softterm2lam))
                  - Ksoftlam * difflam * cos_diff;
      }
    } //ANITA

  }

  return 0;
  
}

int interaction::Perturbed_Soft_Angle_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  m_timer.start();

  SPLIT_VIRIAL_BOUNDARY(_calculate_perturbed_soft_angle_interactions,
			topo, conf, sim);

  m_timer.stop();

  return 0;
}

int interaction::Perturbed_Soft_Angle_Interaction
::init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os,
             bool quiet) 
    {
      if (!quiet)
           os << "Perturbed harmonic soft angle interaction\n";
           
      // add additional angle types with K=0 and the target angle of 
      // the type we are perturbing to or from
      std::vector<topology::perturbed_three_body_term_struct>::iterator bt_it = topo.perturbed_solute().softangles().begin(),
                                 bt_to = topo.perturbed_solute().softangles().end();
      for( ; bt_it != bt_to; ++bt_it){
          if (bt_it->A_type==INT_MAX-1) {
            bt_it->A_type=topo.angle_types_cosharm().size();
            double cost = topo.angle_types_cosharm()[bt_it->B_type].cos0;
            topo.angle_types_cosharm().push_back(interaction::angle_type_struct(0, cost));   
            DEBUG(10, "adding new angle type for soft angle perturbation: " 
                       << bt_it->A_type << " K=" << 0 << " cost=" << cost);         
          } else if (bt_it->B_type==INT_MAX-1) {
            bt_it->B_type=topo.angle_types_cosharm().size();
            double cost = topo.angle_types_cosharm()[bt_it->A_type].cos0;
            topo.angle_types_cosharm().push_back(interaction::angle_type_struct(0, cost)); 
            DEBUG(10, "adding new angle type for soft angle perturbation: " 
                       << bt_it->B_type << " K=" << 0 << " cost=" << cost);      
          }
      }
      return 0;
    }
