/**
 * @file perturbed_nonbonded_term.cc
 * inline methods of Perturbed_Nonbonded_Term
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

/**
 * helper function to initialize the constants.
 */
inline void interaction::Perturbed_Nonbonded_Term
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

  m_crf_cut = (1 - m_crf / 2.0)
    / sim.param().longrange.rf_cutoff;

  // Perturbation
  m_crf_2 = m_crf / 2.0;
  m_cut2 = sim.param().longrange.rf_cutoff * sim.param().longrange.rf_cutoff;

  m_lambda_exp = sim.param().perturbation.lambda_exponent;
  
}

/**
 * helper function to calculate the force and energy for 
 * the reaction field contribution for a given pair
 * using the soft interaction
 */
inline void interaction::Perturbed_Nonbonded_Term
::rf_soft_interaction(math::Vec const &r, double const A_q, double const B_q,
		      double const l, double const alpha_crf,
		      math::Vec & force, double &e_rf, double & de_rf,
		      bool selfterm_correction)
{
  const double dist2 = abs2(r);

  const double A_cut2soft = m_cut2 + alpha_crf * m_B_lambda2;
  const double B_cut2soft = m_cut2 + alpha_crf * m_A_lambda2;

  const double A_cut2soft3 = A_cut2soft * A_cut2soft * A_cut2soft;
  const double B_cut2soft3 = B_cut2soft * B_cut2soft * B_cut2soft;

  const double A_crf_2cut3i = m_crf_2 / sqrt(A_cut2soft3);
  const double B_crf_2cut3i = m_crf_2 / sqrt(B_cut2soft3);

  const double A_crf_cut3i = 2*A_crf_2cut3i;
  const double B_crf_cut3i = 2*B_crf_2cut3i;

  const double A_crf_pert = 3.0*A_crf_2cut3i/A_cut2soft;
  const double B_crf_pert = 3.0*B_crf_2cut3i/B_cut2soft;
 
  force = (m_A_lambda_n * A_q * A_crf_cut3i +
	   m_B_lambda_n * B_q * B_crf_cut3i) * math::four_pi_eps_i *r;

  double const A_e_rf = A_q * (- A_crf_2cut3i * dist2 - m_crf_cut);
  double const B_e_rf = B_q * (- B_crf_2cut3i * dist2 - m_crf_cut);
  
  e_rf = (m_A_lambda_n * A_e_rf + m_B_lambda_n * B_e_rf) * math::four_pi_eps_i;

  if(selfterm_correction)
    e_rf += A_q * math::four_pi_eps_i * m_crf_cut;
  
  DEBUG(11, "A_crf_pert: " << A_crf_pert << " B_crf_pert: " << B_crf_pert);
  DEBUG(11, "m_lambda_exp: " << m_lambda_exp);

  de_rf = ((m_A_lambda_n * A_q * m_B_lambda * A_crf_pert -
	    m_B_lambda_n * B_q * m_A_lambda * B_crf_pert) * dist2 * alpha_crf +
	   (m_B_lambda_n_1 * B_e_rf -
	    m_A_lambda_n_1 * A_e_rf) * m_lambda_exp) * math::four_pi_eps_i;
  
}

inline void interaction::Perturbed_Nonbonded_Term
::rf_scaled_interaction(math::Vec const &r, double const A_q, double const B_q,
			double const l, double const alpha_crf,
			double const A_scale, double const B_scale,
			math::Vec & force, double &e_rf, double & de_rf,
			bool selfterm_correction)
{
  const double dist2 = abs2(r);

  const double A_cut2soft = m_cut2 + alpha_crf * m_B_lambda2;
  const double B_cut2soft = m_cut2 + alpha_crf * m_A_lambda2;

  const double A_cut2soft3 = A_cut2soft*A_cut2soft*A_cut2soft;
  const double B_cut2soft3 = B_cut2soft*B_cut2soft*B_cut2soft;

  const double A_crf_2cut3i = m_crf_2 / sqrt(A_cut2soft3);
  const double B_crf_2cut3i = m_crf_2 / sqrt(B_cut2soft3);

  const double A_crf_cut3i = 2*A_crf_2cut3i;
  const double B_crf_cut3i = 2*B_crf_2cut3i;

  const double A_crf_pert = 3.0*A_crf_2cut3i/A_cut2soft;
  const double B_crf_pert = 3.0*B_crf_2cut3i/B_cut2soft;
 
  const double A_scaled_lambda_n = A_scale * m_A_lambda_n;  
  const double B_scaled_lambda_n = B_scale * m_B_lambda_n;
  
  force = (A_scaled_lambda_n * A_q * A_crf_cut3i +
	   B_scaled_lambda_n * B_q * B_crf_cut3i) * math::four_pi_eps_i *r;

  double const A_e_rf = A_q * (- A_crf_2cut3i * dist2 - m_crf_cut);
  double const B_e_rf = B_q * (- B_crf_2cut3i * dist2 - m_crf_cut);
  
  e_rf = (A_scaled_lambda_n * A_e_rf + B_scaled_lambda_n * B_e_rf) * math::four_pi_eps_i;
  if(selfterm_correction)
    e_rf += A_q * math::four_pi_eps_i * m_crf_cut;
  
  de_rf = ((A_scaled_lambda_n * A_q * m_B_lambda * A_crf_pert -
	    B_scaled_lambda_n * B_q * m_A_lambda * B_crf_pert) * dist2 * alpha_crf +
	   (B_scale * m_B_lambda_n_1 * B_e_rf - 
	    A_scale * m_A_lambda_n_1 * A_e_rf) * m_lambda_exp) * math::four_pi_eps_i;

}

inline void interaction::Perturbed_Nonbonded_Term
::lj_crf_soft_interaction(math::Vec const &r,
			  double const A_c6, double const A_c12,
			  double const B_c6, double const B_c12,
			  double const A_q, double const B_q,
			  double const alpha_lj, double const alpha_crf,
			  double &force1, double &force6, double &force12,
			  double &e_lj, double &e_crf, 
			  double &de_lj, double & de_crf)
{
  double A_c126, B_c126;

  if (A_c6 != 0) A_c126=A_c12/A_c6;
  else A_c126 = 0.0;
  if (B_c6 != 0) B_c126=B_c12/B_c6;
  else B_c126 = 0.0;
  
  const double dist2 = abs2(r);
  assert(dist2 != 0);
  
  const double A_dist2soft = dist2 + alpha_crf*m_B_lambda2;
  const double B_dist2soft = dist2 + alpha_crf*m_A_lambda2;

  const double A_distisoft = 1.0 / sqrt(A_dist2soft);
  const double B_distisoft = 1.0 / sqrt(B_dist2soft);

  const double A_dist3isoft = A_distisoft / A_dist2soft;
  const double B_dist3isoft = B_distisoft / B_dist2soft;
  
  const double dist4 = dist2 * dist2;
  const double dist6 = dist4 * dist2;
  
  const double A_dist6soft = dist6 + alpha_lj*m_B_lambda2*A_c126;
  const double B_dist6soft = dist6 + alpha_lj*m_A_lambda2*B_c126;

  const double A_dist6isoft = 1.0 / A_dist6soft;
  const double B_dist6isoft = 1.0 / B_dist6soft;
  
  const double A_cut2soft = m_cut2 + alpha_crf * m_B_lambda2;
  const double B_cut2soft = m_cut2 + alpha_crf * m_A_lambda2;

  const double A_cut2soft3 = A_cut2soft * A_cut2soft * A_cut2soft;
  const double B_cut2soft3 = B_cut2soft * B_cut2soft * B_cut2soft;

  const double A_crf_2cut3i = m_crf_2 / sqrt(A_cut2soft3);
  const double B_crf_2cut3i = m_crf_2 / sqrt(B_cut2soft3);

  const double A_crf_cut3i = 2 * A_crf_2cut3i;
  const double B_crf_cut3i = 2 * B_crf_2cut3i;

  const double A_crf_pert = 3.0 * A_crf_2cut3i / A_cut2soft;
  const double B_crf_pert = 3.0 * B_crf_2cut3i / B_cut2soft;
  
  
  force1 = (m_A_lambda_n * A_q * (A_distisoft / A_dist2soft + A_crf_cut3i) +
	    m_B_lambda_n * B_q * (B_distisoft / B_dist2soft + B_crf_cut3i)) *
    math::four_pi_eps_i;

  force6 = - 6.0 * (m_A_lambda_n * A_c6 * A_dist6isoft * A_dist6isoft +
		    m_B_lambda_n * B_c6 * B_dist6isoft * B_dist6isoft) * dist4;

  force12 = 12 * (m_A_lambda_n * A_c12 * A_dist6isoft * A_dist6isoft * A_dist6isoft +
		  m_B_lambda_n * B_c12 * B_dist6isoft * B_dist6isoft * B_dist6isoft) * dist4;
  
  const double A_e_lj = (A_c12 * A_dist6isoft - A_c6) * A_dist6isoft;
  const double B_e_lj = (B_c12 * B_dist6isoft - B_c6) * B_dist6isoft;

  const double A_e_crf = A_q * (A_distisoft - A_crf_2cut3i * dist2 - m_crf_cut);
  const double B_e_crf = B_q * (B_distisoft - B_crf_2cut3i * dist2 - m_crf_cut);
  
  DEBUG(11, "just checking\nm_A_lambda " << m_A_lambda
	<< "\nm_B_lambda " << m_B_lambda
	<< "\nm_A_lambda2 " << m_A_lambda2
	<< "\nm_B_lambda2 " << m_B_lambda2
	<< "\nm_A_lambda_n " << m_A_lambda_n
	<< "\nm_B_lambda_n" << m_B_lambda_n
	<< "\nm_A_lambda_n_1 " << m_A_lambda_n_1
	<< "\nm_B_lambda_n_1 " << m_B_lambda_n_1);
  

  e_lj = m_A_lambda_n * A_e_lj + m_B_lambda_n * B_e_lj;
  
  e_crf = (m_A_lambda_n * A_e_crf + m_B_lambda_n * B_e_crf) * math::four_pi_eps_i;
  
  de_lj = -2.0 * alpha_lj * (m_A_lambda_n * m_B_lambda * A_c126 * A_dist6isoft * A_dist6isoft *
			     (2 * A_c12 * A_dist6isoft - A_c6) -
			     m_B_lambda_n * m_A_lambda * B_c126 * B_dist6isoft * B_dist6isoft *
			     (2 * B_c12 * B_dist6isoft - B_c6))
    + m_lambda_exp * (m_B_lambda_n_1 * B_e_lj - m_A_lambda_n_1 * A_e_lj);

  de_crf = -(m_A_lambda_n * A_q * m_B_lambda * (A_dist3isoft - A_crf_pert * dist2) -
	     m_B_lambda_n * B_q * m_A_lambda * (B_dist3isoft - B_crf_pert * dist2)) * 
    math::four_pi_eps_i * alpha_crf
    + m_lambda_exp * (m_B_lambda_n_1 * B_e_crf - m_A_lambda_n_1 * A_e_crf) * math::four_pi_eps_i;
  
}

inline void interaction::Perturbed_Nonbonded_Term
::lj_crf_scaled_interaction(math::Vec const &r,
			    double const A_c6, double const A_c12,
			    double const B_c6, double const B_c12,
			    double const A_q, double const B_q,
			    double const alpha_lj, double const alpha_crf,
			    double const A_scale, double const B_scale,
			    double &force1, double &force6, double &force12,
			    double &e_lj, double &e_crf, 
			    double &de_lj, double & de_crf)
{
  double A_c126, B_c126;

  if (A_c6 != 0) A_c126=A_c12/A_c6;
  else A_c126 = 0.0;
  if (B_c6 != 0) B_c126=B_c12/B_c6;
  else B_c126 = 0.0;
  
  const double dist2 = abs2(r);
  assert(dist2 != 0);
  
  const double A_dist2soft = dist2 + alpha_crf*m_B_lambda2;
  const double B_dist2soft = dist2 + alpha_crf*m_A_lambda2;

  const double A_distisoft = 1.0 / sqrt(A_dist2soft);
  const double B_distisoft = 1.0 / sqrt(B_dist2soft);

  const double A_dist3isoft = A_distisoft / A_dist2soft;
  const double B_dist3isoft = B_distisoft / B_dist2soft;
  
  const double dist4 = dist2 * dist2;
  const double dist6 = dist4 * dist2;
  
  const double A_dist6soft = dist6 + alpha_lj*m_B_lambda2*A_c126;
  const double B_dist6soft = dist6 + alpha_lj*m_A_lambda2*B_c126;

  const double A_dist6isoft = 1.0 / A_dist6soft;
  const double B_dist6isoft = 1.0 / B_dist6soft;
  
  const double A_cut2soft = m_cut2 + alpha_crf * m_B_lambda2;
  const double B_cut2soft = m_cut2 + alpha_crf * m_A_lambda2;

  const double A_cut2soft3 = A_cut2soft * A_cut2soft * A_cut2soft;
  const double B_cut2soft3 = B_cut2soft * B_cut2soft * B_cut2soft;

  const double A_crf_2cut3i = m_crf_2 / sqrt(A_cut2soft3);
  const double B_crf_2cut3i = m_crf_2 / sqrt(B_cut2soft3);

  const double A_crf_cut3i = 2 * A_crf_2cut3i;
  const double B_crf_cut3i = 2 * B_crf_2cut3i;

  const double A_crf_pert = 3.0 * A_crf_2cut3i / A_cut2soft;
  const double B_crf_pert = 3.0 * B_crf_2cut3i / B_cut2soft;
  
  const double A_scaled_lambda_n = A_scale * m_A_lambda_n;
  const double B_scaled_lambda_n = B_scale * m_B_lambda_n;
  
  
  force1 = (A_scaled_lambda_n * A_q * (A_distisoft / A_dist2soft + A_crf_cut3i) +
	    B_scaled_lambda_n * B_q * (B_distisoft / B_dist2soft + B_crf_cut3i)) *
    math::four_pi_eps_i;

  force6 = - 6.0 * (A_scaled_lambda_n * A_c6 * A_dist6isoft * A_dist6isoft +
		    B_scaled_lambda_n * B_c6 * B_dist6isoft * B_dist6isoft) * dist4;

  force12 = 12 * (A_scaled_lambda_n * A_c12 * A_dist6isoft * A_dist6isoft * A_dist6isoft +
		  B_scaled_lambda_n * B_c12 * B_dist6isoft * B_dist6isoft * B_dist6isoft) * dist4;
  
  const double A_e_lj = (A_c12 * A_dist6isoft - A_c6) * A_dist6isoft;
  const double B_e_lj = (B_c12 * B_dist6isoft - B_c6) * B_dist6isoft;

  const double A_e_crf = A_q * (A_distisoft - A_crf_2cut3i * dist2 - m_crf_cut);
  const double B_e_crf = B_q * (B_distisoft - B_crf_2cut3i * dist2 - m_crf_cut);
  

  e_lj = A_scaled_lambda_n * A_e_lj + B_scaled_lambda_n * B_e_lj;
  
  e_crf = (A_scaled_lambda_n * A_e_crf + B_scaled_lambda_n * B_e_crf) * math::four_pi_eps_i;
  
  de_lj = -2.0 * alpha_lj * (A_scaled_lambda_n * m_B_lambda * A_c126 * A_dist6isoft * A_dist6isoft *
			     (2 * A_c12 * A_dist6isoft - A_c6) -
			     B_scaled_lambda_n * m_A_lambda * B_c126 * B_dist6isoft * B_dist6isoft *
			     (2 * B_c12 * B_dist6isoft - B_c6))
    + m_lambda_exp * (B_scale * m_B_lambda_n_1 * B_e_lj - A_scale * m_A_lambda_n_1 * A_e_lj);

  de_crf = (-alpha_crf * (A_scaled_lambda_n * A_q * m_B_lambda * (A_dist3isoft - A_crf_pert * dist2) -
			  B_scaled_lambda_n * B_q * m_A_lambda * (B_dist3isoft - B_crf_pert * dist2))
	    + m_lambda_exp * (B_scale * m_B_lambda_n_1 * B_e_crf - 
			      A_scale * m_A_lambda_n_1 * A_e_crf)) 
    * math::four_pi_eps_i;
  
}

/**
 * Perturbation:
 * lambda value for state A
 */
inline double const interaction::Perturbed_Nonbonded_Term::A_lambda()const
{
  return m_A_lambda;
}

/**
 * Perturbation:
 * lambda value for state B
 */
inline double const interaction::Perturbed_Nonbonded_Term::B_lambda()const
{
  return m_B_lambda;
}

/**
 * Perturbation:
 * lambda value for state A to the power nlam
 */
inline double const interaction::Perturbed_Nonbonded_Term::A_lambda_n()const
{
  return m_A_lambda_n;
}

/**
 * Perturbation:
 * lambda value for state B to the power nlam
 */
inline double const interaction::Perturbed_Nonbonded_Term::B_lambda_n()const
{
  return m_B_lambda_n;
}

/**
 * Perturbation:
 * lambda value for state A to the power nlam-1
 */
inline double const interaction::Perturbed_Nonbonded_Term::A_lambda_n_1()const
{
  return m_A_lambda_n_1;
}

/**
 * Perturbation:
 * lambda value for state B to the power nlam-1
 */
inline double const interaction::Perturbed_Nonbonded_Term::B_lambda_n_1()const
{
  return m_B_lambda_n_1;
}

/**
 * Perturbation:
 * set the lambdas
 */
inline void interaction::Perturbed_Nonbonded_Term::set_lambda(double const l, 
						    int const n)
{
  DEBUG(9, "initializing lambdas: l=" << l
	<< " n=" << n);
  m_A_lambda = 1-l;
  m_A_lambda2 = m_A_lambda * m_A_lambda;

  m_B_lambda = l;
  m_B_lambda2 = m_B_lambda * m_B_lambda;
  
  m_A_lambda_n = pow(m_A_lambda, n);
  m_B_lambda_n = pow(m_B_lambda, n);
  m_A_lambda_n_1 = pow(m_A_lambda, n-1);
  m_B_lambda_n_1 = pow(m_B_lambda, n-1);
  
  DEBUG(11, "\tA:     " << m_A_lambda);
  DEBUG(11, "\tB:     " << m_B_lambda);
  DEBUG(11, "\tA^n:   " << m_A_lambda_n);
  DEBUG(11, "\tB^n:   " << m_B_lambda_n);
  DEBUG(11, "\tA^n-1: " << m_A_lambda_n_1);
  DEBUG(11, "\tB^n-1: " << m_B_lambda_n_1);
}
