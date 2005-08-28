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
  switch(sim.param().force.interaction_function){
  case simulation::lj_crf_func :
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
    break;
    
  case simulation::cgrain_func :
    m_lambda_exp = sim.param().perturbation.lambda_exponent;
    
    // cgrain
    A_cg12= - ((12 + 4) * sim.param().longrange.rf_cutoff) /
      ((pow(sim.param().longrange.rf_cutoff, 12 + 2)) * 
       (sim.param().longrange.rf_cutoff) * 
       (sim.param().longrange.rf_cutoff));
    A_cg6=  - ((6  + 4) * sim.param().longrange.rf_cutoff) /
      ((pow(sim.param().longrange.rf_cutoff, 6  + 2)) * 
       (sim.param().longrange.rf_cutoff) *
       (sim.param().longrange.rf_cutoff));
    A_cg1=  - ((1  + 4) * sim.param().longrange.rf_cutoff) /
      ((pow(sim.param().longrange.rf_cutoff, 1  + 2)) * 
       (sim.param().longrange.rf_cutoff) *
       (sim.param().longrange.rf_cutoff));
    
    B_cg12= ((12 + 3) * sim.param().longrange.rf_cutoff) /
      ((pow(sim.param().longrange.rf_cutoff, 12 + 2)) * 
       (sim.param().longrange.rf_cutoff) * 
       (sim.param().longrange.rf_cutoff) * 
       (sim.param().longrange.rf_cutoff));
    B_cg6=  ((6  + 3) * sim.param().longrange.rf_cutoff) /
      ((pow(sim.param().longrange.rf_cutoff, 6  + 2)) * 
       (sim.param().longrange.rf_cutoff) * 
       (sim.param().longrange.rf_cutoff) * 
       (sim.param().longrange.rf_cutoff));
    B_cg1=  ((1  + 3) * sim.param().longrange.rf_cutoff) /
      ((pow(sim.param().longrange.rf_cutoff, 1  + 2)) * 
       (sim.param().longrange.rf_cutoff) * 
       (sim.param().longrange.rf_cutoff) * 
       (sim.param().longrange.rf_cutoff));
    
    cgrain_eps = sim.param().cgrain.EPS;
    nb_cutoff = sim.param().longrange.rf_cutoff;
    break;
  default:
    io::messages.add("Nonbonded_Innerloop",
		     "interaction function not implemented",
		     io::message::critical);
  }  
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
  
   // substitute A_dist3isoft thing. just like here -- daniel 
  force1 = (m_A_lambda_n * A_q * (A_dist3isoft + A_crf_cut3i) +
	    m_B_lambda_n * B_q * (B_dist3isoft + B_crf_cut3i)) *
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
  
  
  force1 = (A_scaled_lambda_n * A_q * (A_dist3isoft + A_crf_cut3i) +
	    B_scaled_lambda_n * B_q * (B_dist3isoft + B_crf_cut3i)) *
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
 * helper function to calculate the force and energy for
 * a given atom pair in the coarse grain model (perturbed)
 */

inline void interaction::Perturbed_Nonbonded_Term
::cgrain_soft_interaction(math::Vec const &r,
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
  
  const double dist = sqrt(dist2);

  const double A_dist2soft = dist2 + alpha_crf*m_B_lambda2;
  const double B_dist2soft = dist2 + alpha_crf*m_A_lambda2;

  const double A_distisoft = 1.0 / sqrt(A_dist2soft);
  const double B_distisoft = 1.0 / sqrt(B_dist2soft);

  const double A_dist3isoft = A_distisoft / A_dist2soft;
  const double B_dist3isoft = B_distisoft / B_dist2soft;
  
  const double A_dist2soft_cut = nb_cutoff * nb_cutoff + alpha_crf * m_B_lambda2;
  const double B_dist2soft_cut = nb_cutoff * nb_cutoff + alpha_crf * m_A_lambda2;
  
  const double A_distisoft_cut = 1.0 / sqrt(A_dist2soft_cut);
  const double B_distisoft_cut = 1.0 / sqrt(B_dist2soft_cut);

  const double A_dist3isoft_cut = A_distisoft_cut / A_dist2soft_cut;
  const double B_dist3isoft_cut = B_distisoft_cut / B_dist2soft_cut;

  const double dist4 = dist2 * dist2;
  const double dist6 = dist4 * dist2;
  
  const double A_dist6soft = dist6 + alpha_lj*m_B_lambda2*A_c126;
  const double B_dist6soft = dist6 + alpha_lj*m_A_lambda2*B_c126;
  
  const double A_dist6isoft = 1.0 / A_dist6soft;
  const double B_dist6isoft = 1.0 / B_dist6soft;

  // setting the cgrain C constant here
  // for the cgrain, the const C depends on lambda
  
  const double A_distsoft_cut = (pow(nb_cutoff, 6) 
				  + A_c126 * alpha_lj * m_B_lambda2);
  const double B_distsoft_cut = (pow(nb_cutoff, 6) 
				  + B_c126 * alpha_lj * m_A_lambda2);
  
  const double A_dist6isoft_cut = 1.0 / A_distsoft_cut;
  const double B_dist6isoft_cut = 1.0 / B_distsoft_cut;

  // const double A_dist_cut2 = A_distsoft_cut * A_distsoft_cut; 
  // const double B_dist_cut2 = B_distsoft_cut * B_distsoft_cut;
  

  // state A
  A_C_cg12 = A_dist6isoft_cut * A_dist6isoft_cut 
	      - A_cg12  / 3 * nb_cutoff * nb_cutoff * nb_cutoff
	      - B_cg12  / 4 * nb_cutoff * nb_cutoff * nb_cutoff * nb_cutoff;
  A_C_cg6  = A_dist6isoft_cut
	      - A_cg6   / 3 * nb_cutoff * nb_cutoff * nb_cutoff
	      - B_cg6   / 4 * nb_cutoff * nb_cutoff * nb_cutoff * nb_cutoff;
  A_C_cg1 = 1.0 / sqrt((nb_cutoff * nb_cutoff + alpha_crf * m_B_lambda2))
              - A_cg1   / 3 * nb_cutoff * nb_cutoff * nb_cutoff
              - B_cg1   / 4 * nb_cutoff * nb_cutoff * nb_cutoff * nb_cutoff;
  
  // state B
  B_C_cg12 = B_dist6isoft_cut * B_dist6isoft_cut 
	      - A_cg12  / 3 * nb_cutoff * nb_cutoff * nb_cutoff
	      - B_cg12  / 4 * nb_cutoff * nb_cutoff * nb_cutoff * nb_cutoff;
  B_C_cg6  = B_dist6isoft_cut
	      - A_cg6   / 3 * nb_cutoff * nb_cutoff * nb_cutoff
	      - B_cg6   / 4 * nb_cutoff * nb_cutoff * nb_cutoff * nb_cutoff; 
  B_C_cg1 = 1.0 / sqrt((nb_cutoff * nb_cutoff + alpha_crf * m_A_lambda2))
              - A_cg1   / 3 * nb_cutoff * nb_cutoff * nb_cutoff
              - B_cg1   / 4 * nb_cutoff * nb_cutoff * nb_cutoff * nb_cutoff;


  DEBUG(11, "just checking\nm_A_lambda " << m_A_lambda
	<< "\nm_B_lambda " << m_B_lambda
	<< "\nm_A_lambda2 " << m_A_lambda2
	<< "\nm_B_lambda2 " << m_B_lambda2
	<< "\nm_A_lambda_n " << m_A_lambda_n
	<< "\nm_B_lambda_n" << m_B_lambda_n
	<< "\nm_A_lambda_n_1 " << m_A_lambda_n_1
	<< "\nm_B_lambda_n_1 " << m_B_lambda_n_1);
  
  // coulomb 
  force1 = (m_A_lambda_n * A_q * (A_dist3isoft + 
				  A_cg1 * dist + 
				  B_cg1 * dist2) +
	    m_B_lambda_n * B_q * (B_dist3isoft + 
				  A_cg1 * dist + 
				  B_cg1 * dist2)) *
    math::four_pi_eps_i / cgrain_eps;
  
  const double A_e_crf = A_q * (A_distisoft 
				- A_cg1  / 3 * dist2 * dist
				- B_cg1  / 4 * dist2 * dist2
				- A_C_cg1 );
  const double B_e_crf = B_q * (B_distisoft 
				- A_cg1  / 3 * dist2 * dist
				- B_cg1  / 4 * dist2 * dist2
				- B_C_cg1 );
  
  e_crf = (m_A_lambda_n * A_e_crf + m_B_lambda_n * B_e_crf) * math::four_pi_eps_i / cgrain_eps;
  
  de_crf = -(m_A_lambda_n * A_q * m_B_lambda * (A_dist3isoft - A_dist3isoft_cut) -
	     m_B_lambda_n * B_q * m_A_lambda * (B_dist3isoft - B_dist3isoft_cut)) 
    * math::four_pi_eps_i * alpha_crf / cgrain_eps
    + m_lambda_exp * (m_B_lambda_n_1 * B_e_crf - m_A_lambda_n_1 * A_e_crf) 
    * math::four_pi_eps_i / cgrain_eps;

  //LJ
  // after some consideration and an appropriate measure of doubt we
  // change the second 6.0 * into a 1.0 *
  force6 = - 6.0 * (m_A_lambda_n * A_c6 * A_dist6isoft * A_dist6isoft +
		    m_B_lambda_n * B_c6 * B_dist6isoft * B_dist6isoft) * dist4
    -        1.0 * (m_A_lambda_n * A_c6 * (A_cg6  * dist2 + B_cg6  * dist2 * dist) +
		    m_B_lambda_n * B_c6 * (A_cg6  * dist2 + B_cg6  * dist2 * dist)) / dist;
  
  force12 = 12.0 * (m_A_lambda_n * A_c12 * A_dist6isoft * A_dist6isoft * A_dist6isoft +
                    m_B_lambda_n * B_c12 * B_dist6isoft * B_dist6isoft * B_dist6isoft) * dist4 +
    +        1.0 * (m_A_lambda_n * A_c12 * (A_cg12 * dist2 + B_cg12 * dist2 * dist) +
		    m_B_lambda_n * B_c12 * (A_cg12 * dist2 + B_cg12 * dist2 * dist)) / dist;  
  

  const double A_e_lj = A_c12 * (A_dist6isoft * A_dist6isoft 
				 - A_cg12  / 3 * dist2 * dist
				 - B_cg12  / 4 * dist2 * dist2 
				 - A_C_cg12)
    -                    A_c6 * (A_dist6isoft
				 - A_cg6   / 3 * dist2 * dist
				 - B_cg6   / 4 * dist2 * dist2 
				 - A_C_cg6);
  
  const double B_e_lj = B_c12 * (B_dist6isoft * B_dist6isoft 
				 - A_cg12  / 3 * dist2 * dist
				 - B_cg12  / 4 * dist2 * dist2
				 - B_C_cg12)
    -                    B_c6 * (B_dist6isoft
				 - A_cg6   / 3 * dist2 * dist
				 - B_cg6   / 4 * dist2 * dist2
				 - B_C_cg6);
  
  e_lj = m_A_lambda_n * A_e_lj + m_B_lambda_n * B_e_lj; 
  
  de_lj = -2.0 * alpha_lj * (m_A_lambda_n * m_B_lambda * A_c126 * A_dist6isoft * A_dist6isoft *
                             (2 * A_c12 * A_dist6isoft - A_c6) -
                             m_B_lambda_n * m_A_lambda * B_c126 * B_dist6isoft * B_dist6isoft *
                             (2 * B_c12 * B_dist6isoft - B_c6))
    + m_lambda_exp * (m_B_lambda_n_1 * B_e_lj - m_A_lambda_n_1 * A_e_lj);
  
  // std::cout.precision(10); 
  DEBUG(11, "cgrain_pert\nr_ij " << dist
	<< "\nf1 "  << force1
	<< "\ne_c " << e_crf
	<< "\nde_c " << de_crf
	<< "\nf6 " << force6
	<< "\nf12 " << force12
	<< "\nf12_6" << force12+force6
	<< "\ne_lj " << e_lj
	<< "\nde_lj " << de_lj);
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
