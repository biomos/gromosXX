/**
 * @file angle_interaction.tcc
 * template methods of angle_interaction.
 */

/**
 * Constructor.
 */
template<typename t_simulation, typename t_interaction_spec>
inline interaction::angle_interaction<t_simulation, t_interaction_spec>
::angle_interaction()
  : Interaction<t_simulation, t_interaction_spec>("BondAngle")
{
}

/**
 * Destructor.
 */
template<typename t_simulation, typename t_interaction_spec>
inline interaction::angle_interaction<t_simulation, t_interaction_spec>
::~angle_interaction()
{
}

/**
 * calculate angle forces and energies.
 */
template<typename t_simulation, typename t_interaction_spec>
inline void interaction::angle_interaction<t_simulation, t_interaction_spec>
::calculate_interactions(t_simulation &sim)
{
  // loop over the bonds
  std::vector<simulation::Angle>::iterator a_it =
    sim.topology().solute().angles().begin(),
    a_to = sim.topology().solute().angles().end();

  math::VArray &pos   = sim.system().pos();
  math::VArray &force = sim.system().force();
  math::Vec rij, rkj, fi, fj, fk;

  double energy;

  for( ; a_it != a_to; ++a_it){
    
    sim.system().periodicity().nearest_image(pos(a_it->i), pos(a_it->j), rij);
    sim.system().periodicity().nearest_image(pos(a_it->k), pos(a_it->j), rkj);

    double dij = sqrt(dot(rij, rij));
    double dkj = sqrt(dot(rkj, rkj));
    
    assert(dij != 0.0);
    assert(dkj != 0.0);

    double ip = dot(rij, rkj);
    double cost = ip / (dij * dkj);
        
    assert(unsigned(a_it->type) < m_angle_parameter.size());
 
    double K    = m_angle_parameter[a_it->type].K;
    double cos0 = m_angle_parameter[a_it->type].cos0;

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

    if (t_interaction_spec::do_virial == atomic_virial){
      for(int a=0; a<3; ++a)
	for(int b=0; b<3; ++b)
	  sim.system().virial()(a, b) += 
	    rij(a) * fi(b) +
	    rkj(a) * fk(b);

      DEBUG(7, "\tatomic virial done");
    }


    energy = 0.5 * K * (cost - cos0) * (cost - cos0);
    sim.system().energies().angle_energy[sim.topology().
					 atom_energy_group()[a_it->i]]
      += energy;

  }
    
}

/**
 * add angle type.
 */
template<typename t_simulation, typename t_interaction_spec>
inline void interaction::angle_interaction<t_simulation, t_interaction_spec>
::add(angle_type_struct s)
{
  m_angle_parameter.push_back(s);
}

/**
 * add angle type.
 */
template<typename t_simulation, typename t_interaction_spec>
inline void interaction::angle_interaction<t_simulation, t_interaction_spec>
::add(double K, double cos0)
{
  angle_type_struct s;
  s.K = K;
  s.cos0 = cos0;
  add(s);
}

/**
 * access bond parameter.
 */
template<typename t_simulation, typename t_interaction_spec>
inline std::vector<interaction::angle_type_struct> const &
interaction::angle_interaction<t_simulation, t_interaction_spec>
::parameter()const
{
  return m_angle_parameter;
}
