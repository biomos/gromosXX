/**
 * @file nonbonded_interaction.tcc
 * template methods of nonbonded_interaction.
 */

/**
 * calculate nonbonded forces and energies.
 */
template<typename t_simulation>
inline void interaction::nonbonded_interaction<t_simulation>
::calculate_interactions(t_simulation &simu)
{
  double const C12 = 1.0;
  double const C6 =  1.0;
  
  // create a pairlist
  m_pairlist.make_pairlist(simu);
  
  // calculate forces / energies
  math::VArray &pos = simu.system().pos();
  math::VArray &force = simu.system().force();
  
  typename simple_pairlist<t_simulation>::iterator it = m_pairlist.begin();
  for( ; !it.eol(); ++it){
    math::Vec v = pos(it.i()) - pos(it.j());
    double dist2 = dot(v, v);

    assert(dist2 != 0.0);

    double dist6i = 1.0 / (dist2 * dist2 * dist2);
    
    double f = (2 * C12 * dist6i - C6) * 6 * dist6i;
    force(it.i()) += f * v;
    force(it.j()) -= f * v;
  }
  
}
