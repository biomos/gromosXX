/**
 * @file angle_interaction.tcc
 * template methods of angle_interaction.
 */

/**
 * Destructor.
 */
template<typename t_simulation>
inline interaction::angle_interaction<t_simulation>
::~angle_interaction()
{
}

/**
 * calculate angle forces and energies.
 */
template<typename t_simulation>
inline void interaction::angle_interaction<t_simulation>
::calculate_interactions(t_simulation &sim)
{
  // loop over the bonds
  simulation::Angle::iterator a_it =
    sim.topology().solute().angles().begin();

  math::VArray &pos   = sim.system().pos();
  math::VArray &force = sim.system().force();
  math::Vec rij, rkj, fi, fj, fk;

  for( ; !a_it.eol(); ++a_it){
    
    sim.system().periodicity().nearest_image(pos(a_it.i()), pos(a_it.j()), rij);
    sim.system().periodicity().nearest_image(pos(a_it.k()), pos(a_it.j()), rkj);

    double dij = sqrt(dot(rij, rij));
    double dkj = sqrt(dot(rkj, rkj));
    
    assert(dij != 0.0);
    assert(dkj != 0.0);

    double ip = dot(rij, rkj);
    double cost = ip / (dij * dkj);
        
    assert(unsigned(a_it.type()) < m_angle_parameter.size());
 
    double K    = m_angle_parameter[a_it.type()].K;
    double cos0 = m_angle_parameter[a_it.type()].cos0;

    DEBUG(10, "K=" << K << " cos0=" << cos0 << " dij=" << dij << " dkj=" << dkj);

    double ki = -K * (cost - cos0) / dij;
    double kk = -K * (cost - cos0) / dkj;
    
    DEBUG(10, "cost=" << cost << " ki=" << ki << " kk=" << kk);

    fi = ki*(rkj/dkj - rij/dij * cost);
    fk = kk*(rij/dij - rkj/dkj * cost);
    fj = -1.0 * fi - fk;
    
    force(a_it.i()) += fi;
    force(a_it.j()) += fj;
    force(a_it.k()) += fk;
  }
    
}

/**
 * add angle type.
 */
template<typename t_simulation>
inline void interaction::angle_interaction<t_simulation>
::add(angle_type_struct s)
{
  m_angle_parameter.push_back(s);
}

/**
 * add angle type.
 */
template<typename t_simulation>
inline void interaction::angle_interaction<t_simulation>
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
template<typename t_simulation>
inline std::vector<interaction::angle_type_struct> const &
interaction::angle_interaction<t_simulation>
::parameter()const
{
  return m_angle_parameter;
}
