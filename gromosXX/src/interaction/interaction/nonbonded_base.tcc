/**
 * @file nonbonded_base.tcc
 * inline methods of Nonbonded_Base
 */

#undef MODULE
#undef SUBMODULE

#define MODULE interaction
#define SUBMODULE interaction

#include "../../debug.h"

/**
 * Constructor.
 */
inline interaction::Nonbonded_Base::Nonbonded_Base()
{
}

/**
 * add a lj parameter struct.
 */
inline void interaction::Nonbonded_Base
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
inline void interaction::Nonbonded_Base
::coulomb_constant(double const coulomb_constant)
{
  m_coulomb_constant = coulomb_constant;
}

/**
 * get the coulomb constant
 */
inline double interaction::Nonbonded_Base
::coulomb_constant()const
{
  return m_coulomb_constant;
}

/**
 * get the lj parameter for atom types iac_i, iac_j
 */
inline interaction::lj_parameter_struct const &
interaction::Nonbonded_Base
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
inline void interaction::Nonbonded_Base
::resize(size_t i)
{
  m_lj_parameter.resize(i);
  std::vector< std::vector<lj_parameter_struct> >::iterator
    it = m_lj_parameter.begin(),
    to = m_lj_parameter.end();
  
  for(; it!=to; ++it)
    it->resize(i);
}

/**
 * helper function to initialize the constants.
 */
template<typename t_simulation>
inline void interaction::Nonbonded_Base
::initialize(t_simulation const &sim)
{
  // Force
  m_cut3i = 
    1.0 / ( sim.nonbonded().RF_cutoff() 
	    * sim.nonbonded().RF_cutoff() 
	    * sim.nonbonded().RF_cutoff());

  m_crf_cut3i = sim.nonbonded().RF_constant() * m_cut3i;

  // Energy
  m_crf_2cut3i = sim.nonbonded().RF_constant() / 2.0 * m_cut3i;

  m_crf_cut = (1 - sim.nonbonded().RF_constant() / 2.0)
    / sim.nonbonded().RF_cutoff();

}

/**
 * helper function to calculate the force and energy for
 * a given atom pair.
 */
inline void interaction::Nonbonded_Base
::lj_crf_interaction(math::Vec const &r,
		     double const c6, double const c12,
		     double const q,
		     math::Vec &force, double &e_lj, double &e_crf)
{
  assert(dot(r,r) != 0);
  const double dist2 = dot(r, r);
  const double dist2i = 1.0 / dist2;
  const double dist6i = dist2i * dist2i * dist2i;
  const double disti = sqrt(dist2i);
  
  force = ((2 * c12 * dist6i - c6) * 6.0 * dist6i * dist2i + 
    q * coulomb_constant() * (disti * dist2i + m_crf_cut3i)) * r;

  e_lj = (c12 * dist6i - c6) * dist6i;
  e_crf = q * coulomb_constant() * 
    (disti - m_crf_2cut3i * dist2 - m_crf_cut);
  
}

/**
 * helper function to calculate the force and energy for
 * the reaction field contribution for a given pair
 */
inline void interaction::Nonbonded_Base
::rf_interaction(math::Vec const &r,double const q,
		 math::Vec &force, double &e_crf)
{
  const double dist2 = dot(r, r);
  
  force = q * coulomb_constant() *  m_crf_cut3i * r;

  e_crf = q * coulomb_constant() * ( -m_crf_2cut3i * dist2 - m_crf_cut);
  DEBUG(11, "dist2 " << dist2 );
  DEBUG(11, "crf_2cut3i " << m_crf_2cut3i);
  DEBUG(11, "crf_cut " << m_crf_cut);
  DEBUG(11, "q*q   " << q );
  
}

