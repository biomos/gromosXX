/**
 * @file nonbonded_base.tcc
 * inline methods of Nonbonded_Base
 */

#undef MODULE
#undef SUBMODULE

#define MODULE interaction
#define SUBMODULE interaction

#include "../../../debug.h"

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

  // Perturbation
  m_crf_2 = sim.nonbonded().RF_constant() / 2.0;
  m_cut2 = sim.nonbonded().RF_cutoff() * sim.nonbonded().RF_cutoff();
  set_lambda(sim.topology().lambda(), sim.topology().nlam());
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
/**
 * helper function to calculate the force and energy for 
 * the reaction field contribution for a given pair
 * using the soft interaction
 */
inline void interaction::Nonbonded_Base
::rf_soft_interaction(math::Vec const &r, double const q, double const l,
		      double const alpha_crf,
		      math::Vec & force, double &e_rf, double & de_rf)
{
  const double dist2 = dot(r, r);
  const double cut2soft = m_cut2 + alpha_crf*l*l;
  const double cut2soft3 = cut2soft*cut2soft*cut2soft;
  const double crf_2cut3i = m_crf_2 / sqrt(cut2soft3);
  const double crf_cut3i = 2*crf_2cut3i;
  const double crf_pert = 3.0*crf_2cut3i/cut2soft;
 
  force = q * coulomb_constant() * crf_cut3i * r;
  
  e_rf = q * coulomb_constant() * (- crf_2cut3i * dist2 - m_crf_cut);
  de_rf = q*coulomb_constant() * l * alpha_crf * crf_pert*dist2;
  
}


inline void interaction::Nonbonded_Base
::lj_crf_soft_interaction(math::Vec const &r,
			  double const c6, double const c12,
			  double const q, double const l,
			  double const alpha_lj, double const alpha_crf,
			  math::Vec &force, double &e_lj, double &e_crf, 
			  double &de_lj, double & de_crf)
{
  assert(dot(r,r) != 0);

  double c126;
  if (c6 != 0) c126=c12/c6;
  else c126 = 0.0;
  
  const double dist2 = dot(r, r);
  const double dist2soft = dist2 + alpha_crf*l*l;
  const double distisoft = 1.0 / sqrt(dist2soft);
  const double dist3isoft = distisoft / dist2soft;
  
  const double dist6soft = dist2*dist2*dist2 + alpha_lj*l*l*c126;
  const double dist6isoft = 1.0/dist6soft;
  
  const double cut2soft = m_cut2 + alpha_crf*l*l;
  const double cut2soft3 = cut2soft*cut2soft*cut2soft;
  const double crf_2cut3i = m_crf_2 / sqrt(cut2soft3);
  const double crf_cut3i = 2*crf_2cut3i;
  const double crf_pert = 3.0*crf_2cut3i/cut2soft;
  
  
  //const double dist2i = 1.0 / dist2;
  //const double dist6i = dist2i * dist2i * dist2i;
  //const double disti = sqrt(dist2i);
  
  force = ((2 * c12 * dist6isoft - c6) * 
	   6.0 * dist6isoft * dist6isoft * dist2 * dist2 + 
    q * coulomb_constant() * (distisoft / dist2soft + crf_cut3i)) * r;

  e_lj = (c12 * dist6isoft - c6) * dist6isoft;
  e_crf = q * coulomb_constant() * 
    (distisoft - crf_2cut3i * dist2 - m_crf_cut);
  
  de_lj = -2.0 * alpha_lj * l * c126 * dist6isoft * dist6isoft *
    (2 * c12 * dist6isoft - c6);

  de_crf = -q*coulomb_constant() * l * alpha_crf * 
    (dist3isoft - crf_pert*dist2);
  
}

inline double const  interaction::Nonbonded_Base
::crf_2cut3i()const
{
  return m_crf_2cut3i;
}
/**
 * Perturbation:
 * lambda value for state A
 */
inline double const interaction::Nonbonded_Base::A_lambda()const
{
  return m_A_lambda;
}
/**
 * Perturbation:
 * lambda value for state B
 */
inline double const interaction::Nonbonded_Base::B_lambda()const
{
  return m_B_lambda;
}

/**
 * Perturbation:
 * lambda value for state A to the power nlam
 */
inline double const interaction::Nonbonded_Base::A_lambda_n()const
{
  return m_A_lambda_n;
}
/**
 * Perturbation:
 * lambda value for state B to the power nlam
 */
inline double const interaction::Nonbonded_Base::B_lambda_n()const
{
  return m_B_lambda_n;
}
/**
 * Perturbation:
 * lambda value for state A to the power nlam-1
 */
inline double const interaction::Nonbonded_Base::A_lambda_n_1()const
{
  return m_A_lambda_n_1;
}
/**
 * Perturbation:
 * lambda value for state B to the power nlam-1
 */
inline double const interaction::Nonbonded_Base::B_lambda_n_1()const
{
  return m_B_lambda_n_1;
}
/**
 * Perturbation:
 * set the lambdas
 */
inline void interaction::Nonbonded_Base::set_lambda(double const l, 
						    int const n)
{
  DEBUG(5, "initializing lambdas");
  m_A_lambda = 1-l;
  m_B_lambda = l;
  m_A_lambda_n = pow(m_A_lambda, n);
  m_B_lambda_n = pow(m_B_lambda, n);
  m_A_lambda_n_1 = pow(m_A_lambda, n-1);
  m_B_lambda_n_1 = pow(m_B_lambda, n-1);
  DEBUG(5, "\tA:     " << m_A_lambda);
  DEBUG(5, "\tB:     " << m_B_lambda);
  DEBUG(5, "\tA^n:   " << m_A_lambda_n);
  DEBUG(5, "\tB^n:   " << m_B_lambda_n);
  DEBUG(5, "\tA^n-1: " << m_A_lambda_n_1);
  DEBUG(5, "\tB^n-1: " << m_B_lambda_n_1);
}

 
