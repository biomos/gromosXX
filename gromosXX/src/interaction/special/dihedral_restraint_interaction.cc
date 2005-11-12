/**
 * @file dihedral_restraint_interaction.cc
 * template methods of Dihedral_Restraint_Interaction
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

#include <interaction/special/dihedral_restraint_interaction.h>

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
static int _calculate_dihedral_restraint_interactions
(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim)
{
  // loop over the dihedral restraints
  std::vector<topology::dihedral_restraint_struct>::const_iterator 
    it = topo.dihedral_restraints().begin(),
    to = topo.dihedral_restraints().end();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec rij, rkj, rkl, rlj, rmj, rnk, fi, fj, fk, fl;
  double dkj2, dkj, dmj2, dmj, dnk2, dnk, ip, phi;
  double energy, f;

  // math::Vec rij, rkj, rkl, rlj, rim, rln, rmj, rnk, fi, fj, fk, fl;
  // double dkj2, dim, dln, ip;
  // double energy, f;

  math::Periodicity<B> periodicity(conf.current().box);

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
    
    ip = -dot(rij, rnk);
    if(ip < 0) phi *= -1.0;

    /*
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
    */

    DEBUG(9, "uncorrected phi=" << 180.0 * phi / math::Pi);

    while(phi < it->delta)
      phi += 2 * math::Pi;
    while(phi > it->delta + 2 * math::Pi)
      phi -= 2 * math::Pi;

    DEBUG(9, "phi=" << 180 * phi / math::Pi << " phi0=" << 180 * it->phi / math::Pi
	  << " delta=" << 180 * it->delta / math::Pi);

    double delta_phi = phi - it->phi;
    DEBUG(9, "delta_phi=" << 180 * delta_phi / math::Pi);
    
    double phi_lin = sim.param().dihrest.phi_lin;
    double K = sim.param().dihrest.K;
    if (sim.param().dihrest.dihrest == 3)
      K *= it->w0;

    if (phi_lin >= 0.0 && fabs(delta_phi) > phi_lin){
      // LINEAR
      DEBUG(10, "linear");
      double zeta = 1;
      if (delta_phi < 0) zeta = -1;
      
      energy = K * (zeta * delta_phi - 0.5 * phi_lin) * phi_lin;
      f = -K * zeta * phi_lin;
    }
    {
      // HARMONIC
      DEBUG(10, "harmonic");
      energy = 0.5 * K * delta_phi * delta_phi;
      f = -K * delta_phi;
    }

    DEBUG(10, "energy=" << energy << " force=" << f);

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
    
    force(it->i) -= fi;
    force(it->j) -= fj;
    force(it->k) -= fk;
    force(it->l) -= fl;

    /*
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
    */

  }
  
  return 0;
}

int interaction::Dihedral_Restraint_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{

  SPLIT_VIRIAL_BOUNDARY(_calculate_dihedral_restraint_interactions,
			topo, conf, sim);
  
  return 0;
}
