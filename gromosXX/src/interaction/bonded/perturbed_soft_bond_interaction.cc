/**
 * @file perturbed_soft_bond_interaction.cc
 * template methods of Perturbed_Soft_Bond_Interaction
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
#include "perturbed_soft_bond_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE interaction

/**
 * calculate bond forces and energies and lambda derivatives.
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_perturbed_soft_interactions
(  topology::Topology & topo,
   configuration::Configuration & conf,
   simulation::Simulation & sim)
{

  std::vector<interaction::bond_type_struct> const & bondtypes=topo.bond_types_harm();
  DEBUG(7, "perturbed soft bond interaction");
  DEBUG(9, std::setprecision(5));
  
  // loop over the bonds
  std::vector<topology::perturbed_two_body_term_struct>::const_iterator b_it =
    topo.perturbed_solute().softbonds().begin(),
    b_to = topo.perturbed_solute().softbonds().end();
    
  // and the softness parameters
  std::vector<double>::const_iterator alpha_it =
    topo.perturbed_solute().alpha_bond().begin();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec v, f;

  double e = 0.0, diff = 0.0, e_lambda = 0.0;

  math::Periodicity<B> periodicity(conf.current().box);

  for( ; b_it != b_to; ++b_it,++alpha_it){

    // atom i determines the energy group for the output. 
    // we use the same definition for the individual lambdas
    const double lambda = topo.individual_lambda(simulation::bond_lambda)
      [topo.atom_energy_group()[b_it->i]][topo.atom_energy_group()[b_it->i]];
    const double lambda_derivative = topo.individual_lambda_derivative
      (simulation::bond_lambda)
      [topo.atom_energy_group()[b_it->i]][topo.atom_energy_group()[b_it->i]];
       
    DEBUG(7, "soft bond " << b_it->i << "-" << b_it->j
	  << " A-type " << b_it->A_type
	  << " B-type " << b_it->B_type
	  << " lambda " << lambda);

    assert(pos.size() > (b_it->i) && pos.size() > (b_it->j));

    periodicity.nearest_image(pos(b_it->i), pos(b_it->j), v);

    double dist = sqrt(abs2(v));

    DEBUG(7, "dist: " << dist);

    assert(dist != 0.0);
    assert(unsigned(b_it->A_type) < bondtypes.size());
    assert(unsigned(b_it->B_type) < bondtypes.size());

    const double K_A = bondtypes[b_it->A_type].K;
    const double K_B = bondtypes[b_it->B_type].K;
    
    DEBUG(7, "K_A: " << K_A << ", K_B: " << K_B);

    const double r0 = ((1 - lambda) *
		       bondtypes[b_it->A_type].r0 +
		       lambda *
		       bondtypes[b_it->B_type].r0);
    diff = dist - r0;
    const double diff2= diff*diff;
    
    const double soft_A = 1 + *alpha_it * lambda * diff2;
    const double soft_B = 1 + *alpha_it * (1-lambda) * diff2;
    const double soft_A2 = soft_A*soft_A;
    const double soft_B2 = soft_B*soft_B;
    
    DEBUG(9, "r0: " << r0);
    
    //DEBUG(9, "DF " << K * (diff) << "\n" << math::v2s(v));
    
    f = - v * ((1-lambda)*K_A / soft_A2 + lambda*K_B / soft_B2) * diff/dist;

    DEBUG(9, "FORCE: " << math::v2s(f));

    force(b_it->i) += f;
    force(b_it->j) -= f;

    // if (V == math::atomic_virial){
      for(int a=0; a<3; ++a)
	for(int bb=0; bb<3; ++bb)
	  conf.current().virial_tensor(a, bb) += 
	    v(a) * f(bb);

      DEBUG(7, "\tatomic virial done");
      // }

    const double Ksoft = (1-lambda)*K_A / soft_A + lambda*K_B / soft_B;

    e = 0.5 * Ksoft * diff2;

    DEBUG(9, "energy: " << e);

    const double K_diff = bondtypes[b_it->B_type].K -
      bondtypes[b_it->A_type].K;
    DEBUG(9, "K_diff: " << K_diff);
    
    const double b_diff = bondtypes[b_it->B_type].r0 -
      bondtypes[b_it->A_type].r0;
    DEBUG(9, "b_diff: " << b_diff);

    const double softterm1 = 1 + *alpha_it * diff2;
    const double softterm2 = -2 * *alpha_it * lambda * (1-lambda) * diff * b_diff;
    
    e_lambda = lambda_derivative 
       * ( 0.5 * diff2 * ( K_A /  soft_A2  * ((-1) * softterm1 - softterm2) 
                        +  K_B /  soft_B2 * (softterm1 - softterm2))
       - Ksoft * diff * b_diff);

    DEBUG(9, "e_lambda: " << e_lambda);
    
    assert(conf.current().energies.bond_energy.size() >
	   topo.atom_energy_group()[b_it->i]);
    
    conf.current().energies.
      bond_energy[topo.atom_energy_group()
		  [b_it->i]] += e;
    
    assert(conf.current().perturbed_energy_derivatives.bond_energy.size() >
	   topo.atom_energy_group()[b_it->i]);
    
    conf.current().perturbed_energy_derivatives.
      bond_energy[topo.atom_energy_group()
		  [b_it->i]] += e_lambda;

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

        const double b0lam = (1 - lam) * bondtypes[b_it->A_type].r0 
		                         + lam * bondtypes[b_it->B_type].r0;
        double difflam = dist - b0lam;
        double difflam2 = difflam * difflam; 

        const double soft_Alam = 1 + *alpha_it * lam * difflam2;
        const double soft_Blam = 1 + *alpha_it * (1-lam) * difflam2;
        const double soft_Alam2 = soft_Alam*soft_Alam;
        const double soft_Blam2 = soft_Blam*soft_Blam;

        const double Ksoftlam = (1-lam)*K_A / soft_Alam + lam*K_B / soft_Blam;
        const double softterm1lam = 1 + *alpha_it * difflam2;
        const double softterm2lam = -2 * *alpha_it * lam * (1-lam) * difflam * b_diff;

        conf.current().energies.AB_bond[lam_index] += 0.5 * Ksoftlam * difflam2;
        conf.current().perturbed_energy_derivatives.AB_bond[lam_index] += 
                  0.5 * difflam2 * ( K_A /  soft_Alam2  * ((-1) * softterm1lam - softterm2lam) 
                        +  K_B /  soft_Blam2 * (softterm1lam - softterm2lam))
                  - Ksoftlam * difflam * b_diff; 
      }
    } //ANITA
    
  }
  
  return 0;
  
}

int interaction::Perturbed_Soft_Bond_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  m_timer.start();
  
  SPLIT_VIRIAL_BOUNDARY(_calculate_perturbed_soft_interactions,
			topo, conf, sim);
  m_timer.stop();

  return 0;
}

int interaction::Perturbed_Soft_Bond_Interaction
::init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os,
             bool quiet) 
    {
       if (!quiet)
           os << "Perturbed harmonic soft bond interaction\n";
                 
      // add additional bond types with K=0 and the target distance of 
      // the type we are perturbing to or from
      std::vector<topology::perturbed_two_body_term_struct>::iterator bt_it = topo.perturbed_solute().softbonds().begin(),
                                 bt_to = topo.perturbed_solute().softbonds().end();
      for( ; bt_it != bt_to; ++bt_it){
          if (bt_it->A_type==INT_MAX-1) {
            bt_it->A_type=topo.bond_types_harm().size();
            double r = topo.bond_types_harm()[bt_it->B_type].r0;
            topo.bond_types_harm().push_back(interaction::bond_type_struct(0, r)); 
            DEBUG(10, "adding new angle type for soft angle perturbation: " 
                       << bt_it->A_type << " K=" << 0 << " r=" << r);                     
          } else if (bt_it->B_type==INT_MAX-1) {
            bt_it->B_type=topo.bond_types_harm().size();
            double r = topo.bond_types_harm()[bt_it->A_type].r0;
            topo.bond_types_harm().push_back(interaction::bond_type_struct(0, r)); 
            DEBUG(10, "adding new angle type for soft angle perturbation: " 
                       << bt_it->B_type << " K=" << 0 << " r=" << r);              
          }
      }
      return 0;
    }
