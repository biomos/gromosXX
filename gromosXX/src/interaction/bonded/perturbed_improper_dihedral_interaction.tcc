/**
 * @file perturbed_improper_dihedral_interaction.tcc
 * template methods of Perturbed_Improper_Dihedral_Interaction
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
static int _calculate_perturbed_improper_interactions
( topology::Topology & topo,
  configuration::Configuration & conf,
  simulation::Simulation & sim,
  interaction::Improper_Dihedral_Interaction<t_interaction_spec> const & m_interaction)
{
  // this is repeated code from Improper_Dihedral_Interaction !!!

  DEBUG(5, "perturbed improper dihedral interaction");
  DEBUG(7, "using the improper dihedral interaction: " 
	<< m_interaction.name);
  DEBUG(7, std::setprecision(5));
  
  // loop over the angles
  std::vector<topology::perturbed_four_body_term_struct>::const_iterator i_it =
    topo.perturbed_solute().improper_dihedrals().begin(),
    i_to = topo.perturbed_solute().improper_dihedrals().end();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec rij, rkj, rlj, rkl, rmj, rnk, fi, fj, fk, fl;
  double dkj2, dkj, dmj2, dmj, dnk2, dnk, ip, q;
  double energy, e_lambda;

  math::Periodicity<b> periodicity(conf.current().box);

  for( ; i_it != i_to; ++i_it){

    DEBUG(7, "improper dihedral " << i_it->i << "-" << i_it->j << "-" 
	  << i_it->k << "-" << i_it->l
	  << " A-type " << i_it->A_type
	  << " B-type " << i_it->B_type
	  << " lambda " << topo.lambda());

    assert(pos.size() > int(i_it->i) && pos.size() > int(i_it->j) && 
	   pos.size() > int(i_it->k) && pos.size() > int(i_it->l));

    periodicity.nearest_image(pos(i_it->k), pos(i_it->j), rkj);
    periodicity.nearest_image(pos(i_it->i), pos(i_it->j), rij);
    periodicity.nearest_image(pos(i_it->k), pos(i_it->l), rkl);

    rmj = cross(rij, rkj);
    rnk = cross(rkj, rkl);
    
    dkj2 = dot(rkj, rkj);
    dmj2 = dot(rmj, rmj);
    dnk2 = dot(rnk, rnk);
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
    
    assert(unsigned(i_it->A_type) < m_interaction.parameter().size());

    double K    = (1-topo.lambda()) *
      m_interaction.parameter()[i_it->A_type].K +
      topo.lambda() *
      m_interaction.parameter()[i_it->B_type].K;
    double q0 =  (1-topo.lambda()) *
      m_interaction.parameter()[i_it->A_type].q0 +
      topo.lambda() *
      m_interaction.parameter()[i_it->B_type].q0;

    const double K_diff = 
      m_interaction.parameter()[i_it->B_type].K - 
      m_interaction.parameter()[i_it->A_type].K;
    const double q_diff =
      m_interaction.parameter()[i_it->B_type].q0- 
      m_interaction.parameter()[i_it->A_type].q0;
    
    DEBUG(10, "K=" << K << " q0=" << q0 );

    const double ki = -K * (q - q0) * dkj / dmj2;
    const double kl = K * (q - q0) * dkj / dnk2;
    const double kj1 = dot(rij, rkj) / dkj2 - 1.0;
    const double kj2 = dot(rkl, rkj) / dkj2;
    
    fi = ki * rmj;
    fl = kl * rnk;
    fj = kj1 * fi - kj2 * fl;
    fk = -1.0*(fi + fj + fl);
    
    force(i_it->i) += fi;
    force(i_it->j) += fj;
    force(i_it->k) += fk;
    force(i_it->l) += fl;

    if (t_interaction_spec::do_virial == math::atomic_virial){
      periodicity.nearest_image(pos(i_it->l), pos(i_it->j), rlj);

      for(int a=0; a<3; ++a)
	for(int bb=0; bb<3; ++bb)
	  conf.current().virial_tensor(a, bb) += 
	    rij(a) * fi(bb) +
	    rkj(a) * fk(bb) +
	    rlj(a) * fl(bb);

      DEBUG(11, "\tatomic virial done");
    }

    energy = 0.5 * K * (q-q0) * (q-q0);
    
    e_lambda = 0.5 * ( -2.0 * K * q_diff * (q-q0) +
		       K_diff * (q-q0) * (q-q0));
    
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
    
  }

  return 0;
  
}

template<typename t_interaction_spec>
int interaction::Perturbed_Improper_Dihedral_Interaction<t_interaction_spec>
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  switch(conf.boundary_type){
    case math::vacuum :
      return _calculate_perturbed_improper_interactions<math::vacuum, t_interaction_spec>
	(topo, conf, sim, m_interaction);
      break;
    case math::triclinic :
      return _calculate_perturbed_improper_interactions<math::triclinic, t_interaction_spec>
	(topo, conf, sim, m_interaction);
      break;
    case math::rectangular :
      return _calculate_perturbed_improper_interactions<math::rectangular, t_interaction_spec>
	(topo, conf, sim, m_interaction);
      break;
    default:
      throw std::string("Wrong boundary type");
  }
  
}
