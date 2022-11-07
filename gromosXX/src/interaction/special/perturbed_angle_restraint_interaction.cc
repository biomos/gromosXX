/**
 * @file perturbed_angle_restraint_interaction.cc
 * methods of Perturbed_Angle_Restraint_Interaction
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"

#include "../../math/periodicity.h"

// special interactions
#include "../../interaction/interaction_types.h"

#include "../../interaction/special/perturbed_angle_restraint_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

/**
 * calculate angle restraint interactions
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_perturbed_angle_restraint_interactions
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim
 )
{
  DEBUG(5, "perturbed angle restraint interaction");
  // loop over the angle restraints
  std::vector<topology::perturbed_angle_restraint_struct>::const_iterator 
    it = topo.perturbed_angle_restraints().begin(),
    to = topo.perturbed_angle_restraints().end();
    
  std::vector<double>::iterator ene_it = conf.special().pertangleres.energy.begin();
  std::vector<double>::iterator d_it = conf.special().pertangleres.d.begin();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec rij, rkj, fi, fj, fk;
  double energy = 0.0, en_term = 0.0, f = 0.0, energy_derivative = 0.0, dlam_term = 0.0;

  math::Periodicity<B> periodicity(conf.current().box);

  for(; it != to; ++it, ++ene_it, ++d_it){

    // atom i determines the energy group for the output. 
    // we use the same definition for the individual lambdas
    const double l = topo.individual_lambda(simulation::angres_lambda)
      [topo.atom_energy_group()[it->i]]
      [topo.atom_energy_group()[it->i]];
    const double l_deriv = topo.individual_lambda_derivative
      (simulation::angres_lambda)
      [topo.atom_energy_group()[it->i]]
      [topo.atom_energy_group()[it->i]];

    periodicity.nearest_image(pos(it->k), pos(it->j), rkj);
    periodicity.nearest_image(pos(it->i), pos(it->j), rij);

    double dij = sqrt(abs2(rij));
    double dkj = sqrt(abs2(rkj));
    
    DEBUG(15,"dij="<<dij<<" dkj="<<dkj);

    assert(dij != 0.0);
    assert(dkj != 0.0);

    double ip = dot(rij, rkj);
    double cost = ip / (dij * dkj);
    
    if (cost > 1.0) {
      if (cost < 1.0 + math::epsilon) {
        cost = 1.0;
      } else {
        io::messages.add("angle",
                "cos(theta) > 1.0",
                io::message::critical);
      }
    }
    if (cost < -1.0) {
      if (cost > -1.0 - math::epsilon) {
        cost = -1.0;
      } else {
        io::messages.add("angle",
                "cos(theta) < -1.0",
                io::message::critical);
      }
    }
    
    double theta0_A = it->A_theta;
    double cost0_A = cos(it->A_theta);
    double theta0_B = it->B_theta;
    double cost0_B = cos(it->B_theta);
    double cost0 = (1-l) * cost0_A + l * cost0_B;

    double theta  = acos(cost);
    double theta0  = acos(cost0);
    (*d_it) = theta;   
 
    DEBUG(9, "theta=" << 180 * theta / math::Pi << " theta0=" << 180 * acos(cost0) / math::Pi);
    
    double delta_cost = cost - cost0;
    double delta_theta = theta - theta0;
    double K = sim.param().angrest.K;
    double K_A = sim.param().angrest.K;
    double K_B = sim.param().angrest.K;
    
    if (sim.param().angrest.angrest == simulation::angle_restr_inst_weighted){
      K *= (1.0-l) * it->A_w0 + l * it->B_w0;
      K_A *= it->A_w0;
      K_B *= it->B_w0;
    }

    double prefactor = pow(2.0, it->m + it->n) * pow(l, it->n) * pow(1.0-l, it->m);

    // HARMONIC
    en_term = 0.5 * K * delta_cost * delta_cost;
    dlam_term = 0.5 * (K_B - K_A) * delta_cost * delta_cost 
              + K * delta_cost * (cost0_A - cost0_B);
    f = -prefactor * K * delta_cost;
    
    energy = prefactor * en_term;
    (*ene_it) = energy;
    
    conf.current().energies.angrest_energy[topo.atom_energy_group()
					   [it->i]] += energy;
    
    
    fi = f / dij * (rkj/dkj - rij/dij * cost);
    fk = f / dkj * (rij/dij - rkj/dkj * cost);

    fj = -1.0 * fi - fk;
    
    force(it->i) += fi;
    force(it->j) += fj;
    force(it->k) += fk;

    // lambda derivative

    // divide by zero measure
    double dprefndl = 0.0, dprefmdl = 0.0;
    if (it->n==0) dprefndl = 0;
    else dprefndl = it->n * pow(l, it->n-1) * pow(1.0 - l, it->m);
    
    if (it->m == 0) dprefmdl = 0;
    else dprefmdl = it->m * pow(l, it->n) * pow(1.0 - l, it->m-1);

    double dprefdl = pow(2.0, it->m + it->n) * 
      (dprefndl - dprefmdl) * en_term;
    
    double dpotdl = prefactor * dlam_term;

    energy_derivative = l_deriv * (dprefdl + dpotdl);
    
    conf.current().perturbed_energy_derivatives.angrest_energy
      [topo.atom_energy_group()[it->i]] += energy_derivative;

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

        double Klam = (1-lam)*K_A + lam*K_B;
        double cos0lam = (1-lam)*cost0_A + lam*cost0_B;
        double difflam = cost - cos0lam;
        double difflam2 = difflam * difflam;

        double prefactorlam = pow(2.0, it->m + it->n) * pow(lam, it->n) * pow(1.0-lam, it->m);

        double en_termlam = 0.5 * K * delta_cost * delta_cost;
        double dlam_termlam  = 0.5 * (K_B - K_A) * delta_cost * delta_cost 
                  + K * delta_cost * (cost0_A - cost0_B);
        double energylam = prefactorlam * en_termlam;
        
        // lambda derivative
        double dprefndlam = 0.0, dprefmdlam = 0.0;
        // divide by zero measure
        if (it->n==0) dprefndlam = 0;
        else dprefndlam = it->n * pow(lam, it->n-1) * pow(1.0 - lam, it->m);
        
        if (it->m == 0) dprefmdlam = 0;
        else dprefmdlam = it->m * pow(lam, it->n) * pow(1.0 - lam, it->m-1);

        double dprefdlam = pow(2.0, it->m + it->n) * 
          (dprefndlam - dprefmdlam) * en_termlam;
        
        double dpotdlam = prefactorlam * dlam_termlam;
        double energy_derivativlam = (dprefdlam + dpotdlam);

        conf.current().energies.AB_angres[lam_index] += energylam;
        conf.current().perturbed_energy_derivatives.AB_angres[lam_index] += energy_derivativlam;
      }
    } //ANITA
  }
  return 0;
}

int interaction::Perturbed_Angle_Restraint_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{

  SPLIT_VIRIAL_BOUNDARY(_calculate_perturbed_angle_restraint_interactions,
			topo, conf, sim);
  
  return 0;
}
  
 /**
 * initiate angle restraint interactions
 */
template<math::boundary_enum B>
static void _init_angres_data
(topology::Topology & topo,
 configuration::Configuration & conf)
{
   math::Periodicity<B> periodicity(conf.current().box);
   math::VArray &pos   = conf.current().pos;
 
   math::Vec rij, rkj;
    
  for(std::vector<topology::perturbed_angle_restraint_struct>::const_iterator
        it = topo.perturbed_angle_restraints().begin(),
        to = topo.perturbed_angle_restraints().end(); it != to; ++it) {
        
    DEBUG(9, "init: perturbed angle " << it->i << "-" << it->j << "-" << it->k);

    periodicity.nearest_image(pos(it->i), pos(it->j), rij);
    periodicity.nearest_image(pos(it->k), pos(it->j), rkj);
      
    bool warn=false;
    for (int i=0; i<3;  i++) {
        if ((fabs(rij[i]) > conf.current().box(i)[i]*0.45 && abs(rij[i]) < conf.current().box(i)[i]*0.55)
         || (fabs(rkj[i]) > conf.current().box(i)[i]*0.45 && abs(rkj[i]) < conf.current().box(i)[i]*0.55)) {
          warn=true;
        }
    }
    if (warn) {
        std::ostringstream oss;
        oss << "one or more vectors of your angle restraint are\n"
          << "close to half a box length long in your initial structure, \n"
          << "the angle might be calculated from other periodic copies of the atoms than you intended!\n";
          io::messages.add(oss.str(),  "angle_restraint", io::message::warning);
    }
    
    double dij = sqrt(abs2(rij));
    double dkj = sqrt(abs2(rkj));
    
    DEBUG(15,"dij="<<dij<<" dkj="<<dkj);

    assert(dij != 0.0);
    assert(dkj != 0.0);

    double ip = dot(rij, rkj);
    double cost = ip / (dij * dkj);
    
    if (cost > 1.0) {
      if (cost < 1.0 + math::epsilon) {
        cost = 1.0;
      } else {
        io::messages.add("angle",
                "cos(theta) > 1.0",
                io::message::critical);
      }
    }
    if (cost < -1.0) {
      if (cost > -1.0 - math::epsilon) {
        cost = -1.0;
      } else {
        io::messages.add("angle",
                "cos(theta) < -1.0",
                io::message::critical);
      }
    }
    
    double theta  = acos(cost);

    DEBUG(10, "raw theta="<< 180.0 * theta / math::Pi);
             
    conf.special().pertangleres.d.push_back(theta);
    conf.special().pertangleres.energy.push_back(0.0);
  }
}

int interaction::Perturbed_Angle_Restraint_Interaction::init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os,
		     bool quiet) 
{
  
  SPLIT_BOUNDARY(_init_angres_data, topo, conf);
 
  if (!quiet) {
    os << "Perturbed angle restraint interaction";
    os << std::endl;
  }
  return 0;
}
