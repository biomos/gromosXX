/**
 * @file nonbonded_term.cc
 * inline methods of Nonbonded_Term
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

/**
 * helper function to initialize the constants.
 */
inline void interaction::Nonbonded_Term
::init(simulation::Simulation const &sim)
{
  // Force
  m_cut3i = 
    1.0 / ( sim.param().longrange.rf_cutoff
	    * sim.param().longrange.rf_cutoff
	    * sim.param().longrange.rf_cutoff);

  m_crf = 2*(sim.param().longrange.epsilon - sim.param().longrange.rf_epsilon) * 
    (1.0 + sim.param().longrange.rf_kappa * sim.param().longrange.rf_cutoff) -
    sim.param().longrange.rf_epsilon * (sim.param().longrange.rf_kappa  * 
					sim.param().longrange.rf_cutoff *
					sim.param().longrange.rf_kappa  *
					sim.param().longrange.rf_cutoff);

  m_crf /= (sim.param().longrange.epsilon +2* sim.param().longrange.rf_epsilon) *
    (1.0 + sim.param().longrange.rf_kappa * sim.param().longrange.rf_cutoff) +
    sim.param().longrange.rf_epsilon * (sim.param().longrange.rf_kappa  * 
					sim.param().longrange.rf_cutoff *
					sim.param().longrange.rf_kappa  *
					sim.param().longrange.rf_cutoff);
  m_crf_cut3i = m_crf * m_cut3i;

  // Energy
  m_crf_2cut3i = m_crf_cut3i / 2.0;

  m_crf_cut = (1 - m_crf / 2.0)
    / sim.param().longrange.rf_cutoff;
  
}

/**
 * helper function to calculate the force and energy for
 * a given atom pair.
 */
inline void interaction::Nonbonded_Term
::lj_crf_interaction(math::Vec const &r,
		     double const c6, double const c12,
		     double const q,
		     double &force, double &e_lj, double &e_crf)
{
  DEBUG(14, "\t\tnonbonded term");
  
  assert(abs2(r) != 0);
  const double dist2 = abs2(r);
  const double dist2i = 1.0 / dist2;
  const double q_eps = q * math::four_pi_eps_i;
  const double dist6i = dist2i * dist2i * dist2i;
  const double disti = sqrt(dist2i);
  const double c12_dist6i = c12 * dist6i;
  
  e_lj = (c12_dist6i - c6) * dist6i;

  e_crf = q_eps * 
      (disti - m_crf_2cut3i * dist2 - m_crf_cut);

  force = (c12_dist6i + c12_dist6i - c6) * 6.0 * dist6i * dist2i + 
      q_eps * (disti * dist2i + m_crf_cut3i);

  DEBUG(15, "\t\tq=" << q << " 4pie=" << math::four_pi_eps_i 
	<< " crf_cut2i=" << m_crf_cut3i);
  
}

/**
 * helper function to calculate the force and energy for
 * the reaction field contribution for a given pair
 */
inline void interaction::Nonbonded_Term
::rf_interaction(math::Vec const &r,double const q,
		 math::Vec &force, double &e_crf)
{
  const double dist2 = abs2(r);
  
  force = q * math::four_pi_eps_i *  m_crf_cut3i * r;

  e_crf = q * math::four_pi_eps_i * ( -m_crf_2cut3i * dist2 - m_crf_cut);
  DEBUG(11, "dist2 " << dist2 );
  DEBUG(11, "crf_2cut3i " << m_crf_2cut3i);
  DEBUG(11, "crf_cut " << m_crf_cut);
  DEBUG(11, "q*q   " << q );
  
}

inline double const  interaction::Nonbonded_Term
::crf_2cut3i()const
{
  return m_crf_2cut3i;
}

inline void
interaction::Nonbonded_Term::lj_crf_hessian(math::Vec const &r,
				    double const c6, double const c12,
				    double const q,
				    math::Matrix &hess)
{
  const double r2 = math::abs2(r);
  
  const double r4 = r2*r2;
  const double r8 = r4*r4;
  const double r10 = r8*r2;
  const double r14 = r10*r4;
  const double r16 = r8*r8;
    
  // the LENNARD-JONES part
  
  // get the matrix for the first term
  math::dyade(r, r, hess);

  for(int d1=0; d1 < 3; ++d1){
    // first term
    for(int d2=0; d2 < 3; ++d2){
      hess(d1, d2) *= 168.0 * c12 / r16 - 48.0 * c6 / r10;
    }
    // second term
    hess(d1, d1) += 6.0 * c6 / r8 - 12.0 * c12 / r14;
  }

  const double r3 = sqrt(r4 * r2);
  math::Matrix c;
  math::dyade(r, r, c);

  for(int d1=0; d1<3; ++d1){
    for(int d2=0; d2<3; ++d2){
      c(d1, d2) *= 3.0 / r2;
    }
    c(d1, d1) -= 1.0;
  }

  for(int d1=0; d1<3; ++d1){
    for(int d2=0; d2<3; ++d2){
      // first factor
      c(d1, d2) *= q * math::four_pi_eps_i / r3;
    }
    // reaction field term
    c(d1, d1) -= q * math::four_pi_eps_i * m_crf_cut3i;
  }

  for(int d1=0; d1<3; ++d1){
    for(int d2=0; d2<3; ++d2){
      // first factor
      hess(d1, d2) += c(d1, d2);
    }
  }
  
}

