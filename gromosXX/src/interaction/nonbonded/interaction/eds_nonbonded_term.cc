/**
 * @file eds_nonbonded_term.cc
 * inline methods of Eds_Nonbonded_Term
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

#include "simulation/parameter.h"

/**
 * helper function to initialize the constants.
 */
inline void interaction::Eds_Nonbonded_Term
::init(simulation::Simulation const &sim) {
  double crf, crf_2, crf_cut;
  m_crf_cut.clear();
  m_crf.clear();
  m_crf_2.clear();
  
  cgrain_eps.clear();
  switch (sim.param().force.interaction_function) {
    case simulation::lj_crf_func:
      // Force
      if (sim.param().nonbonded.rf_cutoff > 0.0) {
        m_rrf = sim.param().nonbonded.rf_cutoff;
        crf = 2 * (sim.param().nonbonded.epsilon - sim.param().nonbonded.rf_epsilon) *
                (1.0 + sim.param().nonbonded.rf_kappa * m_rrf) -
                sim.param().nonbonded.rf_epsilon * (sim.param().nonbonded.rf_kappa *
                m_rrf * sim.param().nonbonded.rf_kappa * m_rrf);

        crf /= (sim.param().nonbonded.epsilon + 2 * sim.param().nonbonded.rf_epsilon) *
                (1.0 + sim.param().nonbonded.rf_kappa * m_rrf) +
                sim.param().nonbonded.rf_epsilon * (sim.param().nonbonded.rf_kappa *
                m_rrf * sim.param().nonbonded.rf_kappa * m_rrf);

        crf_cut = (1 - crf / 2.0) / m_rrf;
        crf_2 = crf / 2.0;
	
      } else { // infinity case: rf_cutoff == 0
        crf = -1;
        DEBUG(15, "nonbonded term init: m_crf: " << crf);
        m_rrf = 0.0;
        crf_cut = (1 - crf / 2.0) / sim.param().nonbonded.rf_cutoff;
	crf_2 = 0.0;
	
      }
      m_crf_cut.push_back(crf_cut);
      m_crf.push_back(crf);
      m_crf_2.push_back(crf_2);
      
      break;
      
    case simulation::cggromos_func:
      cgrain_eps.push_back(sim.param().cgrain.EPS);
      cgrain_eps.push_back(sim.param().cgrain.EPSM);
      cgrain_eps.push_back(sim.param().nonbonded.epsilon);
      
      for (unsigned int i = 0; i < cgrain_eps.size(); i++) {
        if (sim.param().nonbonded.rf_cutoff > 0.0) {
          m_rrf = sim.param().nonbonded.rf_cutoff;

          crf = 2 * (cgrain_eps[i] - sim.param().nonbonded.rf_epsilon) *
                  (1.0 + sim.param().nonbonded.rf_kappa * m_rrf) -
                  sim.param().nonbonded.rf_epsilon * (sim.param().nonbonded.rf_kappa *
                  m_rrf * sim.param().nonbonded.rf_kappa * m_rrf);
          crf /= (cgrain_eps[i] + 2 * sim.param().nonbonded.rf_epsilon) *
                  (1.0 + sim.param().nonbonded.rf_kappa * m_rrf) +
                  sim.param().nonbonded.rf_epsilon * (sim.param().nonbonded.rf_kappa *
                  m_rrf * sim.param().nonbonded.rf_kappa * m_rrf);

          crf_cut = (1 - crf / 2.0) / m_rrf;
        } else { // infinity case: rf_cutoff == 0
          crf = -1;
          m_rrf = 0.0;
          crf_cut = (1 - crf / 2.0) / sim.param().nonbonded.rf_cutoff;
        }
        m_crf.push_back(crf);
        m_crf_cut.push_back(crf_cut);
      }
      break;
    case simulation::pol_lj_crf_func:
    case simulation::pol_off_lj_crf_func:
    case simulation::cgrain_func:
    default:
      io::messages.add("Nonbonded_Innerloop",
              "interaction function not implemented",
              io::message::critical);
  }
  m_cut2 = sim.param().nonbonded.rf_cutoff * sim.param().nonbonded.rf_cutoff;
  m_lambda_exp = sim.param().perturbation.lambda_exponent;
}

/**
 * Perturbation:
 * set the lambdas
 */
inline void interaction::Eds_Nonbonded_Term::
set_lambda(double const lj_lambda, double const ljs_lambda,
        double const crf_lambda, double const crfs_lambda,
        double const lj_lambda_derivative, double const ljs_lambda_derivative,
        double const crf_lambda_derivative, double const crfs_lambda_derivative,
        int const n) {

  m_A_ljs_lambda = (1 - ljs_lambda) * ljs_lambda_derivative; // multiplication by lambda interaction derivative
  m_A_ljs_lambda2 = (1 - ljs_lambda) * (1 - ljs_lambda);
  m_A_crfs_lambda = (1 - crfs_lambda) * crfs_lambda_derivative; /* multiplication by lambda interaction derivative */
  m_A_crfs_lambda2 = (1 - crfs_lambda) * (1 - crfs_lambda);

  m_B_ljs_lambda = ljs_lambda * ljs_lambda_derivative; /* multiplication by lambda interaction derivative */
  m_B_ljs_lambda2 = ljs_lambda * ljs_lambda;
  m_B_crfs_lambda = crfs_lambda * crfs_lambda_derivative; /* multiplication by lambda interaction derivative */
  m_B_crfs_lambda2 = crfs_lambda * crfs_lambda;

  m_A_lj_lambda_n = pow(1 - lj_lambda, n);
  m_A_crf_lambda_n = pow(1 - crf_lambda, n);
  m_B_lj_lambda_n = pow(lj_lambda, n);
  m_B_crf_lambda_n = pow(crf_lambda, n);
  m_A_lj_lambda_n_1 = pow(1 - lj_lambda, n - 1) * lj_lambda_derivative; /* multiplication by lambda interaction derivative */
  m_A_crf_lambda_n_1 = pow(1 - crf_lambda, n - 1) * crf_lambda_derivative; /* multiplication by lambda interaction derivative */
  m_B_lj_lambda_n_1 = pow(lj_lambda, n - 1) * lj_lambda_derivative; /* multiplication by lambda interaction derivative */
  m_B_crf_lambda_n_1 = pow(crf_lambda, n - 1) * crf_lambda_derivative; /* multiplication by lambda interaction derivative */

  m_n = n;

}



/**
 * helper function to calculate the force and energy for
 * a given atom pair.
 */
inline void interaction::Eds_Nonbonded_Term
::eds_lj_crf_interaction(const double dist2, const double dist6, 
        const double &c6, const double &c12,
        const double &q,
        double const alpha_lj, double const alpha_crf,
        double & force, double & e_nb, unsigned int eps) {
    DEBUG(14, "\t\tnonbonded term");

    double c126 = 0.0;
    if (c6 != 0.0) 
      c126 = c12 / c6;
    
    const double dist6i_alj_c126 = 1.0 / (alpha_lj * c126 + dist6);
    const double q_eps = q * math::four_pi_eps_i;
    const double c12_dist6i = c12 * dist6i_alj_c126;
    const double ac_disti = 1.0 / (sqrt(alpha_crf + dist2));
    const double ac_rrfi = 1.0 / ((alpha_crf + m_rrf*m_rrf)* sqrt(alpha_crf + m_rrf*m_rrf));
    
    double elj = dist6i_alj_c126 * (c12_dist6i - c6);
    double ecrf = q_eps * (ac_disti - (0.5 * m_crf[eps] * ac_rrfi * dist2) - m_crf_cut[eps]);
    e_nb = elj + ecrf;

    double flj = 6.0 * dist2 * dist2 * dist6i_alj_c126 * dist6i_alj_c126 
                 * (2.0 * c12_dist6i - c6);
    double fcrf = q_eps * (ac_disti * ac_disti * ac_disti + m_crf[eps] * ac_rrfi);
    force = flj + fcrf;
    
    DEBUG(15, "nonbonded energy = " << (elj + ecrf) << ", force = " << (flj + fcrf));
}

/**
 * helper function to calculate the force and energy for
 * a given atom pair: one EDS atom and one perturbed atom.
 */
inline void interaction::Eds_Nonbonded_Term
::eds_pert_lj_crf_interaction(const double dist2, const double dist6, 
        const double &A_c6, const double &A_c12,
        const double &B_c6, const double &B_c12,
	const double &A_q,  const double &B_q,
        double const alpha_lj, double const alpha_crf,
	double & force, 
	double & e_lj, double & e_crf, 
        double & de_lj, double & de_crf, unsigned int eps) {
  DEBUG(14, "\t\tnonbonded eds_pert term");
  double A_c126 = 0.0;
  double B_c126 = 0.0;
  
  if (A_c6 != 0.0) 
    A_c126 = A_c12 / A_c6;
  if (B_c6 != 0.0) 
    B_c126 = B_c12 / B_c6;
 
  const double dist4 = dist2 * dist2;

  const double A_dist2soft = dist2 + alpha_crf*m_B_crfs_lambda2;
  const double B_dist2soft = dist2 + alpha_crf*m_A_crfs_lambda2;

  const double A_distisoft = 1.0 / sqrt(A_dist2soft);
  const double B_distisoft = 1.0 / sqrt(B_dist2soft);

  const double A_dist3isoft = A_distisoft / A_dist2soft;
  const double B_dist3isoft = B_distisoft / B_dist2soft;

  const double A_dist6soft = dist6 + alpha_lj * m_B_ljs_lambda2*A_c126;
  const double B_dist6soft = dist6 + alpha_lj * m_A_ljs_lambda2*B_c126;

  const double A_dist6isoft = 1.0 / A_dist6soft;
  const double B_dist6isoft = 1.0 / B_dist6soft;

  const double A_cut2soft = m_cut2 + alpha_crf * m_B_crfs_lambda2;
  const double B_cut2soft = m_cut2 + alpha_crf * m_A_crfs_lambda2;

  const double A_cut2soft3 = A_cut2soft * A_cut2soft * A_cut2soft;
  const double B_cut2soft3 = B_cut2soft * B_cut2soft * B_cut2soft;

  const double A_crf_2cut3i = m_crf_2[eps] / sqrt(A_cut2soft3);
  const double B_crf_2cut3i = m_crf_2[eps] / sqrt(B_cut2soft3);

  const double A_crf_cut3i = 2 * A_crf_2cut3i;
  const double B_crf_cut3i = 2 * B_crf_2cut3i;

  const double A_crf_pert = 3.0 * A_crf_2cut3i / A_cut2soft;
  const double B_crf_pert = 3.0 * B_crf_2cut3i / B_cut2soft;

  DEBUG(14, "first set of parameters computed\nA_dist2soft" << A_dist2soft
	<< "\nA_distisoft " << A_distisoft
	<< "\nA_dist3isoft " << A_dist3isoft
	<< "\nA_dist6soft " << A_dist6soft
	<< "\nA_dist6isoft " << A_dist6isoft
	<< "\nA_cut2soft " << A_cut2soft
	<< "\nA_cut2soft3 " << A_cut2soft3
	<< "\nA_crf_2cut3i " << A_crf_2cut3i
	<< "\nA_crf_cut3i " << A_crf_cut3i
	<< "\nA_crf_pert " << A_crf_pert);
  
  // substitute A_dist3isoft thing. just like here -- daniel
  force = (m_A_crf_lambda_n * A_q * (A_dist3isoft + A_crf_cut3i) +
          m_B_crf_lambda_n * B_q * (B_dist3isoft + B_crf_cut3i)) *
          math::four_pi_eps_i;

  force += -6.0 * (m_A_lj_lambda_n * A_c6 * A_dist6isoft * A_dist6isoft +
          m_B_lj_lambda_n * B_c6 * B_dist6isoft * B_dist6isoft) * dist4;

  force += 12 * (m_A_lj_lambda_n * A_c12 * A_dist6isoft * A_dist6isoft * A_dist6isoft +
          m_B_lj_lambda_n * B_c12 * B_dist6isoft * B_dist6isoft * B_dist6isoft) * dist4;

  DEBUG(14, "forces computed");
  
  DEBUG(11, "just checking\nm_A_ljs_lambda " << m_A_ljs_lambda
	<< "\nm_B_ljs_lambda " << m_B_ljs_lambda
	<< "\nm_A_ljs_lambda2 " << m_A_ljs_lambda2
	<< "\nm_B_ljs_lambda2 " << m_B_ljs_lambda2
	<< "\nm_A_lj_lambda_n " << m_A_lj_lambda_n
	<< "\nm_B_lj_lambda_n" << m_B_lj_lambda_n
	<< "\nm_A_lj_lambda_n_1 " << m_A_lj_lambda_n_1
	<< "\nm_B_lj_lambda_n_1 " << m_B_lj_lambda_n_1
	<< "\nm_A_crfs_lambda " << m_A_crfs_lambda
	<< "\nm_B_crfs_lambda " << m_B_crfs_lambda
	<< "\nm_A_crfs_lambda2 " << m_A_crfs_lambda2
	<< "\nm_B_crfs_lambda2 " << m_B_crfs_lambda2
	<< "\nm_A_crf_lambda_n " << m_A_crf_lambda_n
	<< "\nm_B_crf_lambda_n" << m_B_crf_lambda_n
	<< "\nm_A_crf_lambda_n_1 " << m_A_crf_lambda_n_1
	<< "\nm_B_crf_lambda_n_1 " << m_B_crf_lambda_n_1);
  


  const double A_e_lj = (A_c12 * A_dist6isoft - A_c6) * A_dist6isoft;
  const double B_e_lj = (B_c12 * B_dist6isoft - B_c6) * B_dist6isoft;

  const double A_e_crf = A_q * (A_distisoft - A_crf_2cut3i * dist2 - m_crf_cut[eps]);
  const double B_e_crf = B_q * (B_distisoft - B_crf_2cut3i * dist2 - m_crf_cut[eps]);

  e_lj = m_A_lj_lambda_n * A_e_lj + m_B_lj_lambda_n * B_e_lj;

  e_crf = (m_A_crf_lambda_n * A_e_crf + m_B_crf_lambda_n * B_e_crf) * math::four_pi_eps_i;

  DEBUG(15, "perturbed_nonbonded_term e_lj: " << e_lj << " e_crf: " << e_crf);
  
  de_lj = -2.0 * alpha_lj * (m_A_lj_lambda_n * m_B_ljs_lambda * A_c126 * A_dist6isoft * A_dist6isoft *
          (2 * A_c12 * A_dist6isoft - A_c6) -
          m_B_lj_lambda_n * m_A_ljs_lambda * B_c126 * B_dist6isoft * B_dist6isoft *
          (2 * B_c12 * B_dist6isoft - B_c6))
          + m_lambda_exp * (m_B_lj_lambda_n_1 * B_e_lj - m_A_lj_lambda_n_1 * A_e_lj);

  de_crf = -(m_A_crf_lambda_n * A_q * m_B_crfs_lambda * (A_dist3isoft - A_crf_pert * dist2) -
          m_B_crf_lambda_n * B_q * m_A_crfs_lambda * (B_dist3isoft - B_crf_pert * dist2)) *
          math::four_pi_eps_i * alpha_crf
          + m_lambda_exp * (m_B_crf_lambda_n_1 * B_e_crf - m_A_crf_lambda_n_1 * A_e_crf) *
          math::four_pi_eps_i;

  DEBUG(15, "nonbonded energy = " << (e_lj+e_crf) << ", force = " << (force));
}

/**
 * helper function to calculate the force and energy for
 * the reaction field contribution for a given pair
 */
inline void interaction::Eds_Nonbonded_Term
::eds_rf_interaction(math::Vec const &r,double q, double const alpha_crf,
		     math::Vec &force, double &e_crf, unsigned int eps)
{
  const double dist2 = abs2(r);
  const double ac_rrfi = 1.0 / ((alpha_crf + m_rrf*m_rrf)* sqrt(alpha_crf + m_rrf*m_rrf));
  
  force = q * math::four_pi_eps_i *  m_crf[eps] * ac_rrfi * r;

  e_crf = q * math::four_pi_eps_i * ( -(0.5 * m_crf[eps] * ac_rrfi * dist2) - m_crf_cut[eps]);
  DEBUG(15, "dist2 " << dist2 );
  DEBUG(15, "q*q   " << q );
  DEBUG(15, "rf energy = " << e_crf);
  
}
