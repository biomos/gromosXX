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
 * Constructor.
 */
template<typename t_simulation, typename t_pairlist>
inline interaction::Nonbonded_Interaction<t_simulation, t_pairlist>
::Nonbonded_Interaction()
  : Interaction<t_simulation>("NonBonded")
{
}

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
  DEBUG(15, "Nonbonded_Interaction::add_lj_parameter " 
	<< iac_i << "-" << iac_j);
  
  assert(iac_i < m_lj_parameter.size());
  assert(iac_j < m_lj_parameter.size());
  assert(iac_i < m_lj_parameter[iac_j].size());
  assert(iac_j < m_lj_parameter[iac_i].size());
  
  m_lj_parameter[iac_i][iac_j] = lj;
  m_lj_parameter[iac_j][iac_i] = lj;
}

/** 
 * set the coulomb constant 
 */
template<typename t_simulation, typename t_pairlist>
inline void interaction::Nonbonded_Interaction<t_simulation, t_pairlist>
::coulomb_constant(double const coulomb_constant)
{
  m_coulomb_constant = coulomb_constant;
}

/**
 * get the coulomb constant
 */
template<typename t_simulation, typename t_pairlist>
inline double interaction::Nonbonded_Interaction<t_simulation, t_pairlist>
::coulomb_constant()const
{
  return m_coulomb_constant;
}

/**
 * get the lj parameter for atom types iac_i, iac_j
 */
template<typename t_simulation, typename t_pairlist>
inline interaction::lj_parameter_struct const &
interaction::Nonbonded_Interaction<t_simulation, t_pairlist>
::lj_parameter(size_t iac_i, size_t iac_j)
{
  DEBUG(15, "Nonbonded_Interaction::get_lj_parameter " 
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
  DEBUG(7, "steps " << sim.steps() << " upd " << sim.nonbonded().update());

  if(!(sim.steps() % sim.nonbonded().update())){
    // create a pairlist
    DEBUG(7, "\tupdate the pairlist");
    m_pairlist.update(sim);
  
    // recalc long-range forces
    DEBUG(7, "\tlong range");
    m_longrange_force.resize(sim.system().force().size());
    m_longrange_force = 0.0;
    
    do_interactions(sim, m_pairlist.long_range().begin(),
		    m_pairlist.long_range().end(),
		    longrange);
  }

  // calculate forces / energies
  DEBUG(7, "\tshort range");
  do_interactions(sim, m_pairlist.short_range().begin(),
		  m_pairlist.short_range().end(),
		  shortrange);

  // add long-range force
  sim.system().force() += m_longrange_force;

  // add 1,4 - interactions
  do_14_interactions(sim);

}

/**
 * helper function to calculate forces and energies, 
 * stores them in the arrays pointed to by parameters
 * to make it usable for longrange calculations.
 */
template<typename t_simulation, typename t_pairlist>
inline void interaction::Nonbonded_Interaction<t_simulation, t_pairlist>
::do_interactions(t_simulation &sim, typename t_pairlist::iterator it, 
		  typename t_pairlist::iterator to, 
		  nonbonded_type_enum range)
{
  math::Vec v;
  math::VArray &pos = sim.system().pos();

  math::VArray *force;
  if (range == shortrange) force = &sim.system().force();
  else force = &m_longrange_force;
  

  DEBUG(7, "\tcalculate interactions");
  // Force
  const double cut3i = 
	1.0 / ( sim.nonbonded().RF_cutoff() 
		* sim.nonbonded().RF_cutoff() 
		* sim.nonbonded().RF_cutoff());

  const double crf_di3 = sim.nonbonded().RF_constant() * cut3i;

  // Energy
  const double crf_2 = sim.nonbonded().RF_constant() / 2.0 * cut3i;

  const double crf_di = (1 - sim.nonbonded().RF_constant() / 2.0)
    / sim.nonbonded().RF_cutoff();
  

  for( ; it != to; ++it){
    
    DEBUG(10, "\tpair\t" << it.i() << "\t" << *it);

    sim.system().periodicity().nearest_image(pos(it.i()), pos(*it), v);
    const double dist2 = dot(v, v);
    DEBUG(10, "dist2 = " << dist2);
    assert(dist2 != 0.0);

    const lj_parameter_struct &lj = 
      lj_parameter(sim.topology().iac(it.i()), sim.topology().iac(*it));

    DEBUG(10, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);

    const double dist6i = 1.0 / (dist2 * dist2 * dist2);
    
    const double f_vdw = ((2 * lj.c12 * dist6i - lj.c6) * 6 * dist6i / dist2);

    const double q = sim.topology().charge()(it.i()) 
                   * sim.topology().charge()(*it);
    DEBUG(10, "\tcharge product: " << q);
    DEBUG(10, "coulo=" << coulomb_constant() << " crf_di3=" << crf_di3);
    
    const double f_el = q * coulomb_constant() 
			  * (sqrt(dist6i) + crf_di3 );
    
    DEBUG(10, "\tf_vdw=" << f_vdw << " f_el=" << f_el);

    (*force)(it.i()) += v*(f_vdw + f_el);
    (*force)(*it) -= v*(f_vdw + f_el);

  }
  
}

/**
 * helper function to calculate the forces and energies from the
 * 1,4 interactions.
 */
template<typename t_simulation, typename t_pairlist>
inline void interaction::Nonbonded_Interaction<t_simulation, t_pairlist>
::do_14_interactions(t_simulation &sim)
{
  math::Vec v;
  math::VArray &pos   = sim.system().pos();
  math::VArray &force = sim.system().force();
  
  DEBUG(7, "\tcalculate 1,4-interactions");
  // Force
  const double cut3i = 
	1.0 / ( sim.nonbonded().RF_cutoff() 
		* sim.nonbonded().RF_cutoff() 
		* sim.nonbonded().RF_cutoff());

  const double crf_di3 = sim.nonbonded().RF_constant() * cut3i;

  // Energy
  const double crf_2 = sim.nonbonded().RF_constant() / 2.0 * cut3i;

  const double crf_di = (1 - sim.nonbonded().RF_constant() / 2.0)
    / sim.nonbonded().RF_cutoff();


  std::set<int>::const_iterator it, to;
  
  for(size_t i=0; i<sim.topology().num_solute_atoms(); ++i){
    it = sim.topology().one_four_pair(i).begin();
    to = sim.topology().one_four_pair(i).end();
    
    for( ; it != to; ++it){
      DEBUG(10, "\tpair " << i << " - " << *it);
      
      sim.system().periodicity().nearest_image(pos(i), pos(*it), v);
      const double dist2 = dot(v, v);
      
      DEBUG(10, "\tdist2 = " << dist2);
      assert(dist2 != 0.0);

      const lj_parameter_struct &lj = 
	lj_parameter(sim.topology().iac(i),
		     sim.topology().iac(*it));

      DEBUG(10, "\tlj-parameter cs6=" << lj.cs6 << " cs12=" << lj.cs12);

      const double q = sim.topology().charge()(i) 
	* sim.topology().charge()(*it);

      DEBUG(10, "\tcharge product: " << q);

      const double dist6i = 1.0 / (dist2 * dist2 * dist2);
    
      const double f_vdw = (2 * lj.cs12 * dist6i - lj.cs6)
		  	   * 6 * dist6i / dist2;
      const double f_el = q * coulomb_constant() 
			  * (sqrt(dist6i) + crf_di3 );

      force(i) += v * (f_vdw + f_el);
      force(*it) -= v * (f_vdw + f_el);

    } // loop over 1,4 pairs
  } // loop over solute atoms
}  

/**
 * pairlist accessor
 */
template<typename t_simulation, typename t_pairlist>
t_pairlist & 
interaction::Nonbonded_Interaction<t_simulation, t_pairlist>
::pairlist()
{
  return m_pairlist;
}

  
