/**
 * @file nonbonded_interaction.tcc
 * template methods of nonbonded_interaction.
 */

/**
 * Destructor.
 */
template<typename t_simulation>
inline interaction::nonbonded_interaction<t_simulation>
::~nonbonded_interaction()
{
}

/**
 * add a lj parameter struct.
 */
template<typename t_simulation>
inline void interaction::nonbonded_interaction<t_simulation>
::add_lj_parameter(int i, int j, lj_parameter_struct lj)
{
  assert(i < m_lj_parameter.size());
  assert(j < m_lj_parameter.size());
  assert(i < m_lj_parameter[j].size());
  assert(j < m_lj_parameter[i].size());
  
  m_lj_parameter[i][j] = lj;
  m_lj_parameter[j][i] = lj;
}

/**
 * resize the matrix.
 */
template<typename t_simulation>
inline void interaction::nonbonded_interaction<t_simulation>
::resize(size_t i)
{
  m_lj_parameter.resize(i);
  typename std::vector< std::vector<lj_parameter_struct> >::iterator
    it = m_lj_parameter.begin(),
    to = m_lj_parameter.end();
  
  for(; it!=to; ++it)
    it->resize(i);
}

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
  m_pairlist.update(simu);
  
  // calculate forces / energies
  math::VArray &pos = simu.system().pos();
  math::VArray &force = simu.system().force();
  
  typename simple_pairlist<t_simulation>::iterator it = m_pairlist.begin();
  for( ; it != m_pairlist.end(); ++it){

    math::Vec v = (pos(it.i()) - pos(*it));
    double dist2 = dot(v, v);

    assert(dist2 != 0.0);

    double dist6i = 1.0 / (dist2 * dist2 * dist2);
    
//    math::Vec f = v * ((2 * C12 * dist6i - C6) * 6 * dist6i);
    math::Vec f = 2 * v;

    force(it.i()) += f;
    force(it.j()) -= f;
  }
  
}
