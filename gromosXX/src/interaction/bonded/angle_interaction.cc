/**
 * @file angle_interaction.cc
 * template methods of Angle_Interaction.
 */

/**
 * calculate angle forces and energies.
 */
template<math::boundary_enum b, typename t_interaction_spec>
static int _calculate_angle_interactions(topology::Topology & topo,
				   configuration::Configuration & conf,
				   simulation::Simulation & sim,
				   std::vector<interaction::angle_type_struct> const & param)
{
  // loop over the bonds
  std::vector<topology::three_body_term_struct>::const_iterator a_it =
    topo.solute().angles().begin(),
    a_to = topo.solute().angles().end();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec rij, rkj, fi, fj, fk;

  double energy;

  math::Periodicity<b> periodicity(conf.current().box);
  
  for( ; a_it != a_to; ++a_it){
    
    periodicity.nearest_image(pos(a_it->i), pos(a_it->j), rij);
    periodicity.nearest_image(pos(a_it->k), pos(a_it->j), rkj);

    double dij = sqrt(dot(rij, rij));
    double dkj = sqrt(dot(rkj, rkj));
    
    assert(dij != 0.0);
    assert(dkj != 0.0);

    double ip = dot(rij, rkj);
    double cost = ip / (dij * dkj);
        
    assert(unsigned(a_it->type) < param.size());
 
    double K    = param[a_it->type].K;
    double cos0 = param[a_it->type].cos0;

    DEBUG(10, "K=" << K << " cos0=" << cos0 << " dij=" << dij << " dkj=" << dkj);

    double ki = -K * (cost - cos0) / dij;
    double kk = -K * (cost - cos0) / dkj;
    
    DEBUG(10, "cost=" << cost << " ki=" << ki << " kk=" << kk);

    fi = ki*(rkj/dkj - rij/dij * cost);
    fk = kk*(rij/dij - rkj/dkj * cost);
    fj = -1.0 * fi - fk;
    
    force(a_it->i) += fi;
    force(a_it->j) += fj;
    force(a_it->k) += fk;

    if (t_interaction_spec::do_virial == math::atomic_virial){
      for(int a=0; a<3; ++a)
	for(int bb=0; bb<3; ++bb)
	  conf.current().virial_tensor(a, bb) += 
	    rij(a) * fi(bb) +
	    rkj(a) * fk(bb);

      DEBUG(11, "\tatomic virial done");
    }


    energy = 0.5 * K * (cost - cos0) * (cost - cos0);
    conf.current().energies.angle_energy[topo.atom_energy_group()[a_it->i]]
      += energy;

  }

  return 0;
  
}

template<typename t_interaction_spec>
int interaction::Angle_Interaction<t_interaction_spec>
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  const double start = util::now();
  
  switch(conf.boundary_type){
    case math::vacuum :
      return _calculate_angle_interactions<math::vacuum, t_interaction_spec>
	(topo, conf, sim, m_parameter);
      break;
    case math::triclinic :
      return _calculate_angle_interactions<math::triclinic, t_interaction_spec>
	(topo, conf, sim, m_parameter);
      break;
    case math::rectangular :
      return _calculate_angle_interactions<math::rectangular, t_interaction_spec>
	(topo, conf, sim, m_parameter);
      break;
    default:
      throw std::string("Wrong boundary type");
  }

  m_timing += util::now() - start;
  
}
