/**
 * @file perturbed_dihedral_restraint_interaction.cc
 * methods of Perturbed_Dihedral_Restraint_Interaction
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>
#include <interaction/interaction.h>

#include <math/periodicity.h>

// special interactions
#include <interaction/interaction_types.h>

#include <interaction/special/perturbed_dihedral_restraint_interaction.h>

#include <util/template_split.h>
#include <util/debug.h>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

/**
 * calculate dihedral restraint interactions
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_perturbed_dihedral_restraint_interactions
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim
 )
{
  // loop over the dihedral restraints
  std::vector<topology::perturbed_dihedral_restraint_struct>::const_iterator 
    it = topo.perturbed_dihedral_restraints().begin(),
    to = topo.perturbed_dihedral_restraints().end();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec rij, rkj, rkl, rlj, rmj, rnk, fi, fj, fk, fl;
  double dkj2, dkj, dmj2, dmj, dnk2, dnk, ip, phi;
  double energy, f, dlam;

  /*
    math::VArray &pos   = conf.current().pos;
    math::VArray &force = conf.current().force;
    math::Vec rij, rkj, rkl, rlj, rim, rln, rmj, rnk, fi, fj, fk, fl;
    double dkj2, dim, dln, ip;
    double energy, f, dlam;
  */

  math::Periodicity<B> periodicity(conf.current().box);
  const double l = topo.lambda();

  for(; it != to; ++it){

    periodicity.nearest_image(pos(it->k), pos(it->j), rkj);
    periodicity.nearest_image(pos(it->i), pos(it->j), rij);
    periodicity.nearest_image(pos(it->k), pos(it->l), rkl);
    
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
   
    double acs = ip / (dmj*dnk);
    if (acs > 1.0){
      std::cout << "numerical error?? in : "
		<< it->i << " - " << it->j
		<< " - " << it->k << " - " << it->l
		<< " : " << acs << std::endl;
      
      if (acs < 1.0 + math::epsilon){
	acs = 1.0;
      }
      else{
	io::messages.add("improper dihedral",
			 "acs > 1.0",
			 io::message::critical);
      }
    }
    
    phi  = acos(acs);

    DEBUG(10, "raw phi="<< 180.0 * phi / math::Pi);
    
    ip = dot(rij, rnk);

    if(ip < 0) phi *= -1.0;

    DEBUG(9, "uncorrected phi=" << 180.0 * phi / math::Pi);
    
    while(phi < it->delta)
      phi += 2 * math::Pi;
    while(phi > it->delta + 2 * math::Pi)
      phi -= 2 * math::Pi;

    double phi0_A = it->A_phi;
    double phi0_B = it->B_phi;

    while(phi0_A < it->delta)
      phi0_A += 2 * math::Pi;
    while(phi0_A > it->delta + 2 * math::Pi)
      phi0_A -= 2 * math::Pi;

    while(phi0_B < it->delta)
      phi0_B += 2 * math::Pi;
    while(phi0_B > it->delta + 2 * math::Pi)
      phi0_B -= 2 * math::Pi;
    
    double phi0 = (1-l) * phi0_A + l * phi0_B;
    double delta_phi = phi - phi0;
    double phi_lin = sim.param().dihrest.phi_lin;
    double K = sim.param().dihrest.K;
    double A_K = sim.param().dihrest.K;
    double B_K = sim.param().dihrest.K;
    
    if (sim.param().dihrest.dihrest == 3){
      K *= (1-l) * it->A_w0 + l * it->B_w0;
      A_K *= it->A_w0;
      B_K *= it->B_w0;
    }

    double prefactor = pow(2, it->m + it->n) * pow(l, it->n) * pow(1-l, it->m);

    if (phi_lin >= 0.0 && fabs(delta_phi) > phi_lin){
      // LINEAR
      double zeta = 1;
      if (delta_phi < 0) zeta = -1;
      
      energy = prefactor * K * (zeta * delta_phi - 0.5 * phi_lin) * phi_lin;
      f = - prefactor * K * zeta * phi_lin;
      dlam = 0.5 * phi_lin * ( (B_K - A_K) * (zeta * delta_phi - 0.5 * phi_lin) +
			       K * zeta * (phi0_A - phi0_B));
    }
    else {
      // HARMONIC
      energy = prefactor * 0.5 * K * delta_phi * delta_phi;
      f = -prefactor * K * delta_phi;
      dlam = 0.5 * ( (B_K - A_K) * delta_phi * delta_phi +
		     2 * K * delta_phi * (phi0_A - phi0_B));
    }

    /*
    std::cout << "DIHREST " << it->i << "-" << it->j << "-" << it->k << "-" << it->l
	      << "\t" << 180 * phi0 / math::Pi 
	      << "\t" << 180 * phi / math::Pi
	      << "\tprefactor = " << prefactor
	      << "\tforce = " << f
	      << "\tenergy = " << energy
	      << std::endl;
    */

    conf.current().energies.dihrest_energy[topo.atom_energy_group()
					   [it->i]] += energy;
    
    const double ki = f * dkj / dmj2;
    const double kl = -f * dkj / dnk2;
    const double kj1 = dot(rij, rkj) / dkj2 - 1.0;
    const double kj2 = dot(rkl, rkj) / dkj2;
    
    fi = ki * rmj;
    fl = kl * rnk;
    fj = kj1 * fi - kj2 * fl;
    fk = -1.0*(fi + fj + fl);
    
    force(it->i) += fi;
    force(it->j) += fj;
    force(it->k) += fk;
    force(it->l) += fl;

    // lambda derivative

    // divide by zero measure
    double dprefndl, dprefmdl;
    if (it->n==0) dprefndl = 0;
    else dprefndl = it->n * pow(l, it->n-1) * pow(1 - l, it->m);
    
    if (it->m == 0) dprefmdl = 0;
    else dprefmdl = it->m * pow(l, it->n) * pow(1 - l, it->m-1);

    // CHRIS: This seems to be wrong!
    // energy was already multiplied by prefactor. Here we should only use
    // the restraining part of the energy, before the multiplication with 
    // prefactor
    double dprefdl = pow(2, it->m + it->n) * 
      (dprefndl - dprefmdl) * energy;
    
    double dpotdl = prefactor * dlam;

    conf.current().perturbed_energy_derivatives.dihrest_energy
      [topo.atom_energy_group()[it->i]] += dprefdl + dpotdl;

  }
  
  return 0;
}

int interaction::Perturbed_Dihedral_Restraint_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{

  SPLIT_VIRIAL_BOUNDARY(_calculate_perturbed_dihedral_restraint_interactions,
			topo, conf, sim);
  
  return 0;
}
