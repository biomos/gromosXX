/**
 * @file improper_dihedral_interaction.tcc
 * template methods of Improper_dihedral_interaction.
 */

/**
 * calculate improper dihedral forces and energies.
 */
template<math::boundary_enum b, typename t_interaction_spec>
static int _calculate_improper_interactions(topology::Topology & topo,
					    configuration::Configuration & conf,
					    simulation::Simulation & sim,
					    std::vector<interaction::improper_dihedral_type_struct>
					    const & param)
{
  // loop over the improper dihedrals
  std::vector<topology::four_body_term_struct>::const_iterator i_it =
    topo.solute().improper_dihedrals().begin(),
    i_to = topo.solute().improper_dihedrals().end();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec rij, rkj, rkl, rlj, rmj, rnk, fi, fj, fk, fl;
  double dkj2, dkj, dmj2, dmj, dnk2, dnk, ip, q;
  double energy;
  
  math::Periodicity<b> periodicity(conf.current().box);

  for( ; i_it != i_to; ++i_it){

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
    
    assert(unsigned(i_it->type) < param.size());
 
    const double K  = param[i_it->type].K;
    const double q0 = param[i_it->type].q0;

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
    conf.current().energies.improper_energy[topo.atom_energy_group()[i_it->i]]
      += energy;
    
  }
  return 0;
  
}


template<typename t_interaction_spec>
int interaction::Improper_Dihedral_Interaction<t_interaction_spec>
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  switch(conf.boundary_type){
    case math::vacuum :
      return _calculate_improper_interactions<math::vacuum, t_interaction_spec>
	(topo, conf, sim, m_parameter);
      break;
    case math::triclinic :
      return _calculate_improper_interactions<math::triclinic, t_interaction_spec>
	(topo, conf, sim, m_parameter);
      break;
    case math::rectangular :
      return _calculate_improper_interactions<math::rectangular, t_interaction_spec>
	(topo, conf, sim, m_parameter);
      break;
    default:
      throw std::string("Wrong boundary type");
  }
  
}
