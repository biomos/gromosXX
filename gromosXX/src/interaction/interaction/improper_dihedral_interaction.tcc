/**
 * @file improper_dihedral_interaction.tcc
 * template methods of Improper_dihedral_interaction.
 */

/**
 * Destructor.
 */
template<typename t_simulation>
inline interaction::Improper_dihedral_interaction<t_simulation>
::~Improper_dihedral_interaction()
{
}

/**
 * calculate improper dihedral forces and energies.
 */
template<typename t_simulation>
inline void interaction::Improper_dihedral_interaction<t_simulation>
::calculate_interactions(t_simulation &simu)
{
  // loop over the improper dihedrals
  simulation::Improper_dihedral::iterator i_it =
    simu.topology().solute().improper_dihedrals().begin();

  math::VArray &pos   = simu.system().pos();
  math::VArray &force = simu.system().force();
  math::Vec rij, rkj, rkl, rmj, rnk, fi, fj, fk, fl;
  double dkj2, dkj, dmj2, dmj, dnk2, dnk, ip, q;
  
  for( ; !i_it.eol(); ++i_it){
    rkj = pos(i_it.k()) - pos(i_it.l());
    rij = pos(i_it.i()) - pos(i_it.j());
    rkl = pos(i_it.k()) - pos(i_it.l());
    rmj = cross(rij, rkj);
    rnk = cross(rkj, rkl);
    
    dkj2 = dot(rkj, rkj);
    dmj2 = dot(rmj, rmj);
    dnk2 = dot(rnk, rnk);
    dkj  = sqrt(dkj2);
    dmj  = sqrt(dmj2);
    dnk  = sqrt(dnk2);
    
    assert(dmj != 0.0);
    assert(dnk != 0.0);

    ip = dot(rmj, rnk);
   
    q  = acos(ip / (dmj*dnk));
    ip = dot(rij, rnk);
    if(ip < 0) q *= -1.0;
    
    assert(unsigned(i_it.type()) < m_improper_dihedral_parameter.size());
 
    double K  = m_improper_dihedral_parameter[i_it.type()].K;
    double q0 = m_improper_dihedral_parameter[i_it.type()].q0;

    double ki = -K * (q - q0) * dkj / dmj2;
    double kl = K * (q - q0) * dkj / dnk2;
    double kj1 = dot(rij, rkj) / dkj2 - 1.0;
    double kj2 = dot(rkl, rkj) / dkj2;
    
    fi = ki * rmj;
    fl = kl * rnk;
    fj = kj1 * fi - kj2 * fl;
    fk = -1.0*(fi + fj + fl);
    
    force(i_it.i()) += fi;
    force(i_it.j()) += fj;
    force(i_it.k()) += fk;
    force(i_it.l()) += fl;
  }
}

/**
 * add improper dihedral type.
 */
template<typename t_simulation>
inline void interaction::Improper_dihedral_interaction<t_simulation>
::add(improper_dihedral_type_struct s)
{
  m_improper_dihedral_parameter.push_back(s);
}

/**
 * add improper dihedral type.
 */
template<typename t_simulation>
inline void interaction::Improper_dihedral_interaction<t_simulation>
::add(double K, double q0)
{
  improper_dihedral_type_struct s;
  s.K = K;
  s.q0 = q0;
  add(s);
}

/**
 * access improper dihedral parameter.
 */
template<typename t_simulation>
inline std::vector<interaction::improper_dihedral_type_struct> const &
interaction::Improper_dihedral_interaction<t_simulation>
::parameter()const
{
  return m_improper_dihedral_parameter;
}
