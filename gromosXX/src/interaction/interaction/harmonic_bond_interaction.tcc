/**
 * @file harmonic_bond_interaction.tcc
 * template methods of harmonic_bond_interaction.
 */

/**
 * Destructor.
 */
template<typename t_simulation>
inline interaction::harmonic_bond_interaction<t_simulation>
::~harmonic_bond_interaction()
{
}

/**
 * calculate harmonic bond forces and energies.
 */
template<typename t_simulation>
inline void interaction::harmonic_bond_interaction<t_simulation>
::calculate_interactions(t_simulation &simu)
{
  // loop over the bonds
  simulation::bond::iterator b_it =
    simu.topology().bonds().begin();

  math::VArray &pos   = simu.system().pos();
  math::VArray &force = simu.system().force();
  math::Vec v, f;

  for( ; !b_it.eol(); ++b_it){
    v = pos(b_it.i()) - pos(b_it.j());
    double dist = sqrt(dot(v, v));
    
    assert(dist != 0.0);
    assert(unsigned(b_it.type()) < m_bond_parameter.size());
    
    f = v * (-m_bond_parameter[b_it.type()].K*
	     (dist - m_bond_parameter[b_it.type()].r0) / dist);
    
    force(b_it.i()) += f;
    force(b_it.j()) -= f;
  }
    
}

/**
 * add bond type.
 */
template<typename t_simulation>
inline void interaction::harmonic_bond_interaction<t_simulation>
::add(bond_type_struct s)
{
  m_bond_parameter.push_back(s);
}

/**
 * add bond type.
 */
template<typename t_simulation>
inline void interaction::harmonic_bond_interaction<t_simulation>
::add(double K, double r0)
{
  bond_type_struct s;
  s.K = K;
  s.r0 = r0;
  add(s);
}
