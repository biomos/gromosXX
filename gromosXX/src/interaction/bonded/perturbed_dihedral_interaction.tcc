/**
 * @file perturbed_dihedral_interaction.tcc
 * template methods of Perturbed_Dihedral_Interaction
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE interaction

#include <util/debug.h>

/**
 * calculate angle forces and energies and lambda derivatives.
 */
template<math::boundary_enum b, typename t_interaction_spec>
static int _calculate_perturbed_dihedral_interactions
(  topology::Topology & topo,
   configuration::Configuration & conf,
   simulation::Simulation & sim,
   interaction::Dihedral_Interaction<t_interaction_spec> const & m_interaction)
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

  double dkj2, dim, dln, ip;
  double A_energy, B_energy, energy, e_lambda;

  const double l=topo.lambda();

  math::Periodicity<b> periodicity(conf.current().box);

  for( ; d_it != d_to; ++d_it){

    DEBUG(7, "dihedral " << d_it->i << "-" << d_it->j << "-" 
	  << d_it->k << "-" << d_it->l
	  << " A-type " << d_it->A_type
	  << " B-type " << d_it->B_type
	  << " lambda " << topo.lambda());
    
    assert(pos.size() > int(d_it->i) && pos.size() > int(d_it->j) && 
	   pos.size() > int(d_it->k) && pos.size() > int(d_it->l));

    periodicity.nearest_image(pos(d_it->k), pos(d_it->j), rkj);
    periodicity.nearest_image(pos(d_it->i), pos(d_it->j), rij);
    periodicity.nearest_image(pos(d_it->k), pos(d_it->l), rkl);

    rmj = cross(rij, rkj);
    rnk = cross(rkj, rkl);
    
    dkj2 = dot(rkj, rkj);

    double frim = dot(rij, rkj)/dkj2;
    double frln = dot(rkl, rkj)/dkj2;

    rim = rij - frim * rkj;
    rln = frln * rkj - rkl;
    dim = sqrt(dot(rim, rim));
    dln = sqrt(dot(rln, rln));
    
    ip = dot(rim, rln);

    double cosphi = ip / (dim*dln);
    
    double cosphi2 = cosphi  * cosphi;
    double cosphi3 = cosphi2 * cosphi;
    double cosphi4 = cosphi3 * cosphi;

    assert(unsigned(d_it->A_type) < m_interaction.parameter().size());
    
    double dcosmphi = 0;
    double cosmphi = 0;

    // first state A 
    switch(m_interaction.parameter()[d_it->A_type].m){
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
      default:
        //io::messages.add("dihedral function not implemented for m>6", 
        //               "dihedral_interaction", io::message::error);
        throw std::runtime_error("dihedral type for m=6 not implemented");

    }
    double     K = m_interaction.parameter()[d_it->A_type].K;
    double delta = m_interaction.parameter()[d_it->A_type].pd;

    DEBUG(10, "dihedral K=" << K << " delta=" << delta << " dcos=" << dcosmphi);

    double ki = -K * delta * dcosmphi / dim;
    double kl = -K * delta * dcosmphi / dln;
    double kj1 = frim - 1.0;
    double kj2 = frln;
    
    A_fi = ki * (rln / dln - rim / dim * cosphi);
    A_fl = kl * (rim / dim - rln / dln * cosphi);
    A_fj = kj1 * A_fi - kj2 * A_fl;
    A_fk = -1.0 * (A_fi + A_fj + A_fl);
    

    A_energy = K * (1 + delta * cosmphi);

    // then state B 
    switch(m_interaction.parameter()[d_it->B_type].m){
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
      default:
        //io::messages.add("dihedral function not implemented for m>6", 
        //               "dihedral_interaction", io::message::error);
        throw std::runtime_error("dihedral type for m=6 not implemented");

    }
        K = m_interaction.parameter()[d_it->B_type].K;
    delta = m_interaction.parameter()[d_it->B_type].pd;

    DEBUG(10, "dihedral K=" << K << " delta=" << delta << " dcos=" << dcosmphi);

    ki = -K * delta * dcosmphi / dim;
    kl = -K * delta * dcosmphi / dln;
    kj1 = frim - 1.0;
    kj2 = frln;
    
    B_fi = ki * (rln / dln - rim / dim * cosphi);
    B_fl = kl * (rim / dim - rln / dln * cosphi);
    B_fj = kj1 * B_fi - kj2 * B_fl;
    B_fk = -1.0 * (B_fi + B_fj + B_fl);
    
    B_energy = K * (1 + delta * cosmphi);

    // now combine
    force(d_it->i) += (1.0-l) * A_fi + l * B_fi;
    force(d_it->j) += (1.0-l) * A_fj + l * B_fj;
    force(d_it->k) += (1.0-l) * A_fk + l * B_fk;
    force(d_it->l) += (1.0-l) * A_fl + l * B_fl;

    if (t_interaction_spec::do_virial == math::atomic_virial){
      periodicity.nearest_image(pos(d_it->l), pos(d_it->j), rlj);

      for(int a=0; a<3; ++a)
	for(int bb=0; bb<3; ++bb)
	  conf.current().virial_tensor(a, bb) += 
	    rij(a) * ((1.0-l) * A_fi(bb) + l * B_fi(bb)) +
	    rkj(a) * ((1.0-l) * A_fk(bb) + l * B_fk(bb)) +
	    rlj(a) * ((1.0-l) * A_fl(bb) + l * B_fl(bb));

      DEBUG(11, "\tatomic virial done");
    }

    energy = (1.0-l) * A_energy + l * B_energy;
    e_lambda = (B_energy - A_energy);

    assert(conf.current().energies.dihedral_energy.size() >
	   topo.atom_energy_group()[d_it->i]);
    conf.current().energies.dihedral_energy
      [topo.atom_energy_group()[d_it->i]] += energy;

   assert(conf.current().perturbed_energy_derivatives.dihedral_energy.size() >
	   topo.atom_energy_group()[d_it->i]);
   conf.current().perturbed_energy_derivatives.dihedral_energy
     [topo.atom_energy_group()[d_it->i]] += e_lambda;   
  }

  return 0;

}

template<typename t_interaction_spec>
int interaction::Perturbed_Dihedral_Interaction<t_interaction_spec>
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation & sim)
{
  switch(conf.boundary_type){
    case math::vacuum :
      return _calculate_perturbed_dihedral_interactions<math::vacuum, t_interaction_spec>
	(topo, conf, sim, m_interaction);
      break;
    case math::triclinic :
      return _calculate_perturbed_dihedral_interactions<math::triclinic, t_interaction_spec>
	(topo, conf, sim, m_interaction);
      break;
    case math::rectangular :
      return _calculate_perturbed_dihedral_interactions<math::rectangular, t_interaction_spec>
	(topo, conf, sim, m_interaction);
      break;
    default:
      throw std::string("Wrong boundary type");
  }
  
}
