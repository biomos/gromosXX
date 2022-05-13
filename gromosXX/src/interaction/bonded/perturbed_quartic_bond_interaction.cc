/**
 * @file perturbed_quartic_bond_interaction.cc
 * template methods of Perturbed_Quartic_Bond_Interaction
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
#include "quartic_bond_interaction.h"
#include "perturbed_quartic_bond_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE bonded

/**
 * calculate quartic bond forces and energies and lambda derivatives.
 */
template<math::boundary_enum B, math::virial_enum V>
int _calculate_perturbed_qbond_interactions
( topology::Topology & topo,
  configuration::Configuration & conf,
  simulation::Simulation & sim,
  interaction::Quartic_Bond_Interaction const & m_interaction)
{
  std::vector<interaction::bond_type_struct> const & bondtypes=topo.bond_types_quart();
  DEBUG(7, "perturbed quartic bond interaction");
  DEBUG(8, "using the bond interaction: " << m_interaction.name);
  DEBUG(8, std::setprecision(5));
  
  // loop over the bonds
  std::vector<topology::perturbed_two_body_term_struct>::iterator b_it =
    topo.perturbed_solute().bonds().begin(),
    b_to = topo.perturbed_solute().bonds().end();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec v, f;

  double e = 0.0, e_lambda = 0.0;

  math::Periodicity<B> periodicity(conf.current().box);

  for( ; b_it != b_to; ++b_it){

    // atom i determines the energy group for the output. 
    // we use the same definition for the individual lambdas
    const double lambda = topo.individual_lambda(simulation::bond_lambda)
      [topo.atom_energy_group()[b_it->i]][topo.atom_energy_group()[b_it->i]];
    const double lambda_derivative = topo.individual_lambda_derivative
      (simulation::bond_lambda)
      [topo.atom_energy_group()[b_it->i]][topo.atom_energy_group()[b_it->i]];
    
    DEBUG(7, "bond " << b_it->i << "-" << b_it->j
	  << " A-type " << b_it->A_type
	  << " B-type " << b_it->B_type 
	  << " lambda " << lambda);

    assert(pos.size() > (b_it->i) && pos.size() > (b_it->j));
    periodicity.nearest_image(pos(b_it->i), pos(b_it->j), v);

    double dist2 = abs2(v);

    DEBUG(7, "dist2: " << dist2);

    assert(unsigned(b_it->A_type) < bondtypes.size());
    assert(unsigned(b_it->B_type) < bondtypes.size());

    const double K = (1 - lambda) * bondtypes[b_it->A_type].K +
      lambda * bondtypes[b_it->B_type].K;
    
    DEBUG(7, "K: " << K);

    const double r0 = ((1 - lambda) *
		       bondtypes[b_it->A_type].r0 +
		       lambda *
		       bondtypes[b_it->B_type].r0);

    const double r02 = r0 * r0;

    DEBUG(7, "r02: " << r02);
    
    DEBUG(7, "DF " << K * (dist2 - r02) << "\n" << math::v2s(v));
    
    f = v * (-K) * (dist2 - r02);

    DEBUG(7, "FORCE: " << math::v2s(f));
    
    force(b_it->i) += f;
    force(b_it->j) -= f;

    // if (V == math::atomic_virial){
      for(int a=0; a<3; ++a)
	for(int c=0; c<3; ++c)
	  conf.current().virial_tensor(a, c) += 
	    v(a) * f(c);

      DEBUG(7, "\tatomic virial done");
      // }
  
    e = 0.25 * K * (dist2 -r02) * (dist2 - r02);

    DEBUG(7, "energy: " << e);

    const double K_diff = bondtypes[b_it->B_type].K -
      bondtypes[b_it->A_type].K;
    DEBUG(7, "K_diff: " << K_diff);
    
    const double b_diff = bondtypes[b_it->B_type].r0 -
      bondtypes[b_it->A_type].r0;

    DEBUG(7, "b_diff: " << b_diff);
    
    const double b_mix = bondtypes[b_it->A_type].r0 +
      lambda * b_diff;
    DEBUG(7, "b_mix: " << b_mix);
    
    e_lambda = 0.25 * lambda_derivative * 
      ( -4 * (bondtypes[b_it->A_type].K +
	      lambda * K_diff) *
	b_diff * b_mix *
	(dist2 - b_mix * b_mix) +
	K_diff * 
	(dist2 - b_mix * b_mix) * (dist2 - b_mix * b_mix));
    DEBUG(7, "e_lambda: " << e_lambda);
    
    assert(conf.current().energies.bond_energy.size() >
	   topo.atom_energy_group()[b_it->i]);
    
    conf.current().energies.
      bond_energy[topo.atom_energy_group()[b_it->i]] += e;
    
    assert(conf.current().perturbed_energy_derivatives.bond_energy.size() >
	   topo.atom_energy_group()[b_it->i]);
    
    conf.current().perturbed_energy_derivatives.
      bond_energy[topo.atom_energy_group()[b_it->i]] += e_lambda;

    // ANITA
    if (sim.param().precalclam.nr_lambdas &&
        ((sim.steps() % sim.param().write.free_energy) == 0)){
      double KA = bondtypes[b_it->A_type].K;
      double KB = bondtypes[b_it->B_type].K;
      double b0A = bondtypes[b_it->A_type].r0;
      double b0B = bondtypes[b_it->B_type].r0;

      double lambda_step = (sim.param().precalclam.max_lam -
                            sim.param().precalclam.min_lam) /
                            (sim.param().precalclam.nr_lambdas-1);

      //loop over nr_lambdas
      for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){

        // determine current lambda for this index
        double lam=(lam_index * lambda_step) + sim.param().precalclam.min_lam;

        double b0lam = (1-lam)*b0A + lam*b0B;
        double b0lam2 = b0lam * b0lam;
        double Klam = (1-lam)*KA + lam*KB; 
        double difflam = dist2 -b0lam2;
        double difflam2 = difflam*difflam;

        conf.current().energies.AB_bond[lam_index] += 0.25 * Klam * difflam2;
        conf.current().perturbed_energy_derivatives.AB_bond[lam_index] += 
                  0.25 * K_diff * difflam2 - Klam * b_diff * b0lam * difflam;
      }
    } //ANITA
  }

  return 0;
    
}

int interaction::Perturbed_Quartic_Bond_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  m_timer.start();
  
  SPLIT_VIRIAL_BOUNDARY(_calculate_perturbed_qbond_interactions,
			topo, conf, sim, m_interaction);

  m_timer.stop();
  return 0;
  
}
