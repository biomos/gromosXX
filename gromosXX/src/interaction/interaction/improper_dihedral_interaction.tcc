/**
 * @file improper_dihedral_interaction.tcc
 * template methods of Improper_dihedral_interaction.
 */

/**
 * Constructor.
 */
template<typename t_simulation>
inline interaction::Improper_dihedral_interaction<t_simulation>
::Improper_dihedral_interaction()
  : Interaction<t_simulation>("ImproperDihedral")
{
}

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
::calculate_interactions(t_simulation &sim)
{
  // loop over the improper dihedrals
  simulation::Improper_dihedral::iterator i_it =
    sim.topology().solute().improper_dihedrals().begin();

  math::VArray &pos   = sim.system().pos();
  math::VArray &force = sim.system().force();
  math::Vec rij, rkj, rkl, rmj, rnk, fi, fj, fk, fl;
  double dkj2, dkj, dmj2, dmj, dnk2, dnk, ip, q;
  double energy;
  
  for( ; !i_it.eol(); ++i_it){
    sim.system().periodicity().nearest_image(pos(i_it.k()), pos(i_it.j()), rkj);
    sim.system().periodicity().nearest_image(pos(i_it.i()), pos(i_it.j()), rij);
    sim.system().periodicity().nearest_image(pos(i_it.k()), pos(i_it.l()), rkl);

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
    
    assert(unsigned(i_it.type()) < m_improper_dihedral_parameter.size());
 
    const double K  = m_improper_dihedral_parameter[i_it.type()].K;
    const double q0 = m_improper_dihedral_parameter[i_it.type()].q0;

    const double ki = -K * (q - q0) * dkj / dmj2;
    const double kl = K * (q - q0) * dkj / dnk2;
    const double kj1 = dot(rij, rkj) / dkj2 - 1.0;
    const double kj2 = dot(rkl, rkj) / dkj2;
    
    fi = ki * rmj;
    fl = kl * rnk;
    fj = kj1 * fi - kj2 * fl;
    fk = -1.0*(fi + fj + fl);
    
    force(i_it.i()) += fi;
    force(i_it.j()) += fj;
    force(i_it.k()) += fk;
    force(i_it.l()) += fl;
    
    energy = 0.5 * K * (q-q0) * (q-q0);
    sim.system().energies().improper_energy[sim.topology().atom_energy_group()[i_it.i()]] += energy;

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
