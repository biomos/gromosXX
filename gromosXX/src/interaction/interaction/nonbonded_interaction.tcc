/**
 * @file nonbonded_interaction.tcc
 * template methods of Nonbonded_Interaction.
 */

#undef MODULE
#undef SUBMODULE

#define MODULE interaction
#define SUBMODULE interaction

#include "../../debug.h"

/**
 * Destructor.
 */
template<typename t_simulation, typename t_pairlist>
inline interaction::Nonbonded_Interaction<t_simulation, t_pairlist>
::~Nonbonded_Interaction()
{
  DEBUG(4, "Nonbonded_Interaction::destructor");
}

/**
 * add a lj parameter struct.
 */
template<typename t_simulation, typename t_pairlist>
inline void interaction::Nonbonded_Interaction<t_simulation, t_pairlist>
::add_lj_parameter(size_t iac_i, size_t iac_j, lj_parameter_struct lj)
{
  DEBUG(4, "Nonbonded_Interaction::add_lj_parameter " 
	<< iac_i << "-" << iac_j);
  
  assert(iac_i < m_lj_parameter.size());
  assert(iac_j < m_lj_parameter.size());
  assert(iac_i < m_lj_parameter[iac_j].size());
  assert(iac_j < m_lj_parameter[iac_i].size());
  
  m_lj_parameter[iac_i][iac_j] = lj;
  m_lj_parameter[iac_j][iac_i] = lj;
}

/**
 * get the lj parameter for atom types iac_i, iac_j
 */
template<typename t_simulation, typename t_pairlist>
inline interaction::lj_parameter_struct const &
interaction::Nonbonded_Interaction<t_simulation, t_pairlist>
::lj_parameter(size_t iac_i, size_t iac_j)
{
  DEBUG(4, "Nonbonded_Interaction::get_lj_parameter " 
	<< iac_i << "-" << iac_j);

  assert(iac_i < m_lj_parameter.size());
  assert(iac_j < m_lj_parameter[iac_i].size());
  
  return m_lj_parameter[iac_i][iac_j];
}

/**
 * resize the matrix.
 */
template<typename t_simulation, typename t_pairlist>
inline void interaction::Nonbonded_Interaction<t_simulation, t_pairlist>
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
template<typename t_simulation, typename t_pairlist>
inline void interaction::Nonbonded_Interaction<t_simulation, t_pairlist>
::calculate_interactions(t_simulation &sim)
{
  DEBUG(4, "Nonbonded_Interaction::calculate_interactions");

  // need to update pairlist?
  if(!(sim.steps() % sim.nonbonded_update())){
    // create a pairlist
    DEBUG(7, "\tupdate the pairlist");
    m_pairlist.update(sim);
  
    // recalc long-range forces
    DEBUG(7, "\tlong range");
    m_longrange_force.resize(sim.system().force().size());
    m_longrange_force = 0.0;
    
    do_interactions(sim, m_pairlist.long_range().begin(),
		    m_pairlist.long_range().end(),
		    m_longrange_force);
  }

  // calculate forces / energies
  DEBUG(7, "\tshort range");
  do_interactions(sim, m_pairlist.short_range().begin(),
		  m_pairlist.short_range().end(),
		  sim.system().force());

  // add long-range force
  sim.system().force() += m_longrange_force;

}

/**
 * helper function to calculate forces and energies, 
 * stores them in the arrays pointed to by parameters
 * to make it usable for longrange calculations.
 */
template<typename t_simulation, typename t_pairlist>
inline void interaction::Nonbonded_Interaction<t_simulation, t_pairlist>
::do_interactions(t_simulation &sim, typename t_pairlist::iterator it, 
		  typename t_pairlist::iterator to, math::VArray &force)
{
  math::Vec v;
  math::VArray &pos = sim.system().pos();

  DEBUG(7, "\tcalculate interactions");

  for( ; it != to; ++it){
    
    DEBUG(10, "\tpair\t" << it.i() << "\t" << *it);

    sim.system().periodicity().nearest_image(pos(it.i()), pos(*it), v);
    const double dist2 = dot(v, v);

    DEBUG(10, "\tdist2 = " << dist2);
    assert(dist2 != 0.0);

    const lj_parameter_struct &lj = 
      lj_parameter(sim.topology().iac(it.i()), sim.topology().iac(*it));

    DEBUG(10, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);

    const double dist6i = 1.0 / (dist2 * dist2 * dist2);
    
    math::Vec f = v * ((2 * lj.c12 * dist6i - lj.c6) * 6 * dist6i);

    force(it.i()) += f;
    force(it.j()) -= f;
  }
  
}

  
