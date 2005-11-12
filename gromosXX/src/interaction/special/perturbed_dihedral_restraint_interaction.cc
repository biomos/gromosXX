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
  math::Vec rij, rkj, rkl, rlj, rim, rln, rmj, rnk, fi, fj, fk, fl;
  double dkj2, dim, dln, ip;
  double energy, f, dlam;

  math::Periodicity<B> periodicity(conf.current().box);
  double l = topo.lambda();

  for(; it != to; ++it){

    periodicity.nearest_image(pos(it->i), pos(it->j), rij);
    periodicity.nearest_image(pos(it->k), pos(it->j), rkj);
    periodicity.nearest_image(pos(it->k), pos(it->l), rkl);
    
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
    double phi = acos(cosphi);
    
    double sgn = dot(rij, rnk);
    if(sgn < 0) phi *= -1.0;

    while(phi < it->delta)
      phi += 2 * math::Pi;
    while(phi > it->delta + 2 * math::Pi)
      phi -= 2 * math::Pi;

    double phi0 = (1-l) * it->A_phi + l * it->B_phi;
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
			       K * zeta * (it->A_phi - it->B_phi));
    }
    {
      // HARMONIC
      energy = prefactor * 0.5 * K * delta_phi * delta_phi;
      f = -prefactor * K * delta_phi;
      dlam = 0.5 * ( (B_K - A_K) * delta_phi * delta_phi +
		     2 * K * delta_phi * (it->A_phi - it->B_phi));
    }

    conf.current().energies.dihrest_energy[topo.atom_energy_group()
					   [it->i]] += energy;
    
    double ki = f / dim;
    double kl = f / dln;
    double kj1 = frim - 1.0;
    double kj2 = frln;
    
    fi = ki * (rln / dln - rim / dim * cosphi);
    fl = kl * (rim / dim - rln / dln * cosphi);
    fj = kj1 * fi - kj2 * fl;
    fk = -1.0 * (fi + fj + fl);
    
    force(it->i) += fi;
    force(it->j) += fj;
    force(it->k) += fk;
    force(it->l) += fl;

    // lambda derivative
    double dprefdl = pow(2, it->m + it->n) * 
      (it->n * pow(l, it->n-1) * pow(1 - l, it->m) - it->m * pow(l, it->n) * pow(1 - l, it->m-1)) * energy;
    
    double dpotdl = prefactor * dlam;

    conf.current().perturbed_energy_derivatives.dihrest_energy[topo.atom_energy_group()[it->i]] += dprefdl + dpotdl;

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
