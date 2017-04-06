/**
 * @file perturbed_new_bond_interaction.cc
 * template methods of Perturbed_New_Bond_Interaction
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
#include "perturbed_new_bond_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE bonded

/**
 * calculate new bond forces and energies and lambda derivatives.
 */
template<math::boundary_enum B, math::virial_enum V>
int _calculate_perturbed_new_bond_interactions
( topology::Topology & topo,
  configuration::Configuration & conf,
  simulation::Simulation & sim,
  interaction::Quartic_Bond_Interaction const & m_interaction)
{

  DEBUG(7, "perturbed new bond interaction");
  DEBUG(8, "using the bond interaction: " << m_interaction.name);
  DEBUG(8, std::setprecision(5));
  
  // loop over the bonds
  std::vector<topology::perturbed_two_body_term_struct>::iterator b_it =
    topo.perturbed_solute().bonds().begin(),
    b_to = topo.perturbed_solute().bonds().end();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec v, f;

  double e, e_lambda;

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

    double dist = abs(v);
    double dist2 = abs2(v);
    int powlamb = sim.param().force.powlamb;

    DEBUG(7, "dist2: " << dist2);
    DEBUG(7, "powlamb: " << powlamb);

    assert(unsigned(b_it->A_type) < m_interaction.parameter().size());
    assert(unsigned(b_it->B_type) < m_interaction.parameter().size());

    double KA = m_interaction.parameter()[b_it->A_type].K;
    double KB = m_interaction.parameter()[b_it->B_type].K;
    double r0A= m_interaction.parameter()[b_it->A_type].r0;
    double r0A2 = r0A * r0A;
    double rdiffA = dist2 - r0A2;
    double rdiffA2 = rdiffA * rdiffA;
    double r0B= m_interaction.parameter()[b_it->B_type].r0;
    double r0B2 = r0B * r0B;
    double rdiffB = dist2 - r0B2;
    double rdiffB2 = rdiffB * rdiffB;
    double powlamX = pow(lambda,powlamb);
    double pow1lamX = pow((1-lambda),powlamb);
    double VA,VB;

    if(r0A==r0B){
      VA = 0.25 * KA * rdiffA2;
      VB = 0.25 * KB * rdiffB2;
      e = (1-lambda) * VA + lambda * VB;
      f = (-1)*((1-lambda) * (v*KA*rdiffA) + lambda*(v*KB*rdiffB));
      e_lambda = -VA + VB;
    }
    else if(((r0A > r0B) && (dist < r0B)) || ((r0B > r0A) && (dist > r0B))){
      VA = 0.25 * pow1lamX * KA * rdiffA2;
      VB = 0.25 * KB * rdiffB2;
      e = (1-lambda) * VA + lambda * VB;
      f = (-1)*((1-lambda) * (v*pow1lamX*KA*rdiffA) + lambda*(v*KB*rdiffB));
      e_lambda = -VA + VB - 3*pow1lamX*KA*rdiffA2;
    }
    else if(((r0A > r0B) && (dist > r0A)) || ((r0B > r0A) && (dist < r0A))){
      VA = 0.25 * KA * rdiffA2;
      VB = 0.25 * powlamX * KB * rdiffB2;
      e = (1-lambda) * VA + lambda * VB;
      f = (-1)*((1-lambda) * (v*KA*rdiffA) + lambda*(v*powlamX*KB*rdiffB));
      e_lambda = -VA + VB + 3*powlamX*KB*rdiffB2;
    }
    else{
      VA = 0.25 * pow1lamX * KA * rdiffA2;
      VB = 0.25 * powlamX * KB * rdiffB2;
      e = (1-lambda) * VA + lambda * VB;
      f = (-1)*((1-lambda) * (v*pow1lamX*KA*rdiffA) + lambda*(v*powlamX*KB*rdiffB));
      e_lambda = -VA + VB - 3*pow1lamX*KA*rdiffA2 + 3*powlamX*KB*rdiffB2;
    }
    
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

      double lambda_step = (sim.param().precalclam.max_lam -
                            sim.param().precalclam.min_lam) /
                            (sim.param().precalclam.nr_lambdas-1);

      if(r0A==r0B){
        VA = 0.25 * KA * rdiffA2;
        VB = 0.25 * KB * rdiffB2;
        e_lambda = -VA + VB;
        // for each lambda:
        for (int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){
          double lam=(lam_index * lambda_step) + sim.param().precalclam.min_lam;
          conf.current().energies.AB_bond[lam_index] += (1-lam) * VA + lam * VB;
          conf.current().perturbed_energy_derivatives.AB_bond[lam_index] += e_lambda;
        }
      }
      else if(((r0A > r0B) && (dist < r0B)) || ((r0B > r0A) && (dist > r0B))){
        // for each lambda:
        for (int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){
          double lam=(lam_index * lambda_step) + sim.param().precalclam.min_lam;
          double powlamX = pow(lam,powlamb);
          double pow1lamX = pow((1-lam),powlamb);
          VA = 0.25 * pow1lamX * KA * rdiffA2;
          VB = 0.25 * KB * rdiffB2;
          conf.current().energies.AB_bond[lam_index] += (1-lam) * VA + lam * VB;
          conf.current().perturbed_energy_derivatives.AB_bond[lam_index] += 
              -VA + VB - 3*pow1lamX*KA*rdiffA2;
        }
      }
      else if(((r0A > r0B) && (dist > r0A)) || ((r0B > r0A) && (dist < r0A))){
        // for each lambda:
        for (int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){
          double lam=(lam_index * lambda_step) + sim.param().precalclam.min_lam;
          double powlamX = pow(lam,powlamb);
          double pow1lamX = pow((1-lam),powlamb);
          VA = 0.25 * KA * rdiffA2;
          VB = 0.25 * powlamX * KB * rdiffB2;
          conf.current().energies.AB_bond[lam_index] += (1-lam) * VA + lam * VB;
          conf.current().perturbed_energy_derivatives.AB_bond[lam_index] += 
               -VA + VB + 3*powlamX*KB*rdiffB2;
        }
      }
      else{
        // for each lambda:
        for (int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){
          double lam=(lam_index * lambda_step) + sim.param().precalclam.min_lam;
          double powlamX = pow(lam,powlamb);
          double pow1lamX = pow((1-lam),powlamb);
          VA = 0.25 * pow1lamX * KA * rdiffA2;
          VB = 0.25 * powlamX * KB * rdiffB2;
          conf.current().energies.AB_bond[lam_index] += (1-lam) * VA + lam * VB;
          conf.current().perturbed_energy_derivatives.AB_bond[lam_index] +=
             -VA + VB - 3*pow1lamX*KA*rdiffA2 + 3*powlamX*KB*rdiffB2;
        }
      }


    } //ANITA
  }

  return 0;
    
}

int interaction::Perturbed_New_Bond_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  m_timer.start();
  
  SPLIT_VIRIAL_BOUNDARY(_calculate_perturbed_new_bond_interactions,
			topo, conf, sim, m_interaction);

  m_timer.stop();
  return 0;
  
}
