/**
 * @file perturbed_nonbonded_term.cc
 * inline methods of Perturbed_Nonbonded_Term
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

#include "simulation/parameter.h"

/**
 * helper function to initialize the constants.
 */
inline void interaction::Perturbed_Nonbonded_Term
::init(simulation::Simulation const &sim) {
    double cut3i = 0.0, crf = 0.0, crf_cut3i = 0.0, crf_cut = 0.0, crf_2 = 0.0;
    m_cut3i.clear();
    m_crf.clear();
    m_crf_cut3i.clear();
    m_crf_cut.clear();
    m_crf_2.clear();
    cgrain_eps.clear();
    switch (sim.param().force.interaction_function) {
        case simulation::lj_crf_func:
        case simulation::lj_crf_2d_func:
        case simulation::pol_lj_crf_func:
        case simulation::pol_off_lj_crf_func:

            // Force
            if (sim.param().nonbonded.rf_cutoff > 0.0) {
                if (sim.param().nonbonded.rf_epsilon == 0.0) {
                    cut3i = 1.0 / (sim.param().nonbonded.rf_cutoff
                            * sim.param().nonbonded.rf_cutoff
                            * sim.param().nonbonded.rf_cutoff);
                    DEBUG(15, "nonbonded term init: m_cut3i: " << cut3i);

                    crf = -1;
                    DEBUG(15, "nonbonded term init: m_crf: " << crf);
                    crf_cut3i = crf * cut3i;

                    // Energy
                    crf_cut = (1 - crf / 2.0)
                            / sim.param().nonbonded.rf_cutoff;
                } else {
                    cut3i = 1.0 / (sim.param().nonbonded.rf_cutoff
                            * sim.param().nonbonded.rf_cutoff
                            * sim.param().nonbonded.rf_cutoff);

                    crf = 2 * (sim.param().nonbonded.epsilon - sim.param().nonbonded.rf_epsilon) *
                            (1.0 + sim.param().nonbonded.rf_kappa * sim.param().nonbonded.rf_cutoff) -
                            sim.param().nonbonded.rf_epsilon * (sim.param().nonbonded.rf_kappa *
                            sim.param().nonbonded.rf_cutoff *
                            sim.param().nonbonded.rf_kappa *
                            sim.param().nonbonded.rf_cutoff);

                    crf /= (sim.param().nonbonded.epsilon + 2 * sim.param().nonbonded.rf_epsilon) *
                            (1.0 + sim.param().nonbonded.rf_kappa * sim.param().nonbonded.rf_cutoff) +
                            sim.param().nonbonded.rf_epsilon * (sim.param().nonbonded.rf_kappa *
                            sim.param().nonbonded.rf_cutoff *
                            sim.param().nonbonded.rf_kappa *
                            sim.param().nonbonded.rf_cutoff);
                    crf_cut3i = crf * cut3i;

                    crf_cut = (1 - crf / 2.0) / sim.param().nonbonded.rf_cutoff;

                    // Perturbation
                }
                crf_2 = crf / 2.0;
                m_cut3i.push_back(cut3i);
                m_crf.push_back(crf);
                m_crf_cut3i.push_back(crf_cut3i);
                m_crf_cut.push_back(crf_cut);
                m_crf_2.push_back(crf_2);
            } else {
                m_cut3i.push_back(0.0);
                m_crf.push_back(0.0);
                m_crf_cut3i.push_back(0.0);
                m_crf_cut.push_back(0.0);
                m_crf_2.push_back(0.0);
            }
            break;
        case simulation::cgrain_func:
            A_cg12 = -(12 * (12 + 4)) / (pow(sim.param().nonbonded.rf_cutoff, 12 + 3));
            A_cg6 = -(6 * (6 + 4)) / (pow(sim.param().nonbonded.rf_cutoff, 6 + 3));
            A_cg1 = -(1 * (1 + 4)) / (pow(sim.param().nonbonded.rf_cutoff, 1 + 3));

            B_cg12 = (12 * (12 + 3)) / (pow(sim.param().nonbonded.rf_cutoff, 12 + 4));
            B_cg6 = (6 * (6 + 3)) / (pow(sim.param().nonbonded.rf_cutoff, 6 + 4));
            B_cg1 = (1 * (1 + 3)) / (pow(sim.param().nonbonded.rf_cutoff, 1 + 4));

            cgrain_eps.push_back(sim.param().cgrain.EPS);
            nb_cutoff = sim.param().nonbonded.rf_cutoff;
            break;
        case simulation::cggromos_func:
            cgrain_eps.push_back(sim.param().cgrain.EPS); // CG-CG
            cgrain_eps.push_back(sim.param().cgrain.EPSM); // FG-CG
            cgrain_eps.push_back(sim.param().nonbonded.epsilon); // FG-FG

            if (sim.param().nonbonded.rf_cutoff > 0.0) {

                cut3i = 1.0 / (sim.param().nonbonded.rf_cutoff
                        * sim.param().nonbonded.rf_cutoff
                        * sim.param().nonbonded.rf_cutoff);
                for (unsigned int i = 0; i < cgrain_eps.size(); ++i) {
                        crf = 2 * (cgrain_eps[i] - sim.param().nonbonded.rf_epsilon) *
                                (1.0 + sim.param().nonbonded.rf_kappa * sim.param().nonbonded.rf_cutoff) -
                                sim.param().nonbonded.rf_epsilon * (sim.param().nonbonded.rf_kappa *
                                sim.param().nonbonded.rf_cutoff *
                                sim.param().nonbonded.rf_kappa *
                                sim.param().nonbonded.rf_cutoff);

                        crf /= (cgrain_eps[i] + 2 * sim.param().nonbonded.rf_epsilon) *
                                (1.0 + sim.param().nonbonded.rf_kappa * sim.param().nonbonded.rf_cutoff) +
                                sim.param().nonbonded.rf_epsilon * (sim.param().nonbonded.rf_kappa *
                                sim.param().nonbonded.rf_cutoff *
                                sim.param().nonbonded.rf_kappa *
                                sim.param().nonbonded.rf_cutoff);
                        crf_cut3i = crf * cut3i;

                        crf_cut = (1 - crf / 2.0) / sim.param().nonbonded.rf_cutoff;
                    // Perturbation
                        crf_2 = crf / 2.0;
                m_cut3i.push_back(0.0);
                m_crf.push_back(0.0);
                m_crf_cut3i.push_back(0.0);
                m_crf_cut.push_back(0.0);
                m_crf_2.push_back(crf_2);
                }
            } else {
                crf_2 = crf / 2.0;
                m_cut3i.push_back(0.0);
                m_crf.push_back(0.0);
                m_crf_cut3i.push_back(0.0);
                m_crf_cut.push_back(0.0);
                m_crf_2.push_back(crf_2);
            }
            
            break;
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
 * lambda value for lj interaction for state A
 */
inline double const & interaction::Perturbed_Nonbonded_Term::A_lj_lambda()const {
  return m_A_lj_lambda;
}

/**
 * Perturbation:
 * lambda value for crf interaction for state A
 */
inline double const & interaction::Perturbed_Nonbonded_Term::A_crf_lambda()const {
  return m_A_crf_lambda;
}

/**
 * Perturbation:
 * lambda value for lj interaction for state B
 */
inline double const & interaction::Perturbed_Nonbonded_Term::B_lj_lambda()const {
  return m_B_lj_lambda;
}

/**
 * Perturbation:
 * lambda value for crf interaction for state B
 */
inline double const & interaction::Perturbed_Nonbonded_Term::B_crf_lambda()const {
  return m_B_crf_lambda;
}

/**
 * Perturbation:
 * lambda value for lj interaction for state A to the power nlam
 */
inline double const & interaction::Perturbed_Nonbonded_Term::A_lj_lambda_n()const {
  return m_A_lj_lambda_n;
}

/**
 * Perturbation:
 * lambda value for crf interaction for state A to the power nlam
 */
inline double const & interaction::Perturbed_Nonbonded_Term::A_crf_lambda_n()const {
  return m_A_crf_lambda_n;
}

/**
 * Perturbation:
 * lambda value for lj interaction for state B to the power nlam
 */
inline double const & interaction::Perturbed_Nonbonded_Term::B_lj_lambda_n()const {
  return m_B_lj_lambda_n;
}

/**
 * Perturbation:
 * lambda value for crf interaction for state B to the power nlam
 */
inline double const & interaction::Perturbed_Nonbonded_Term::B_crf_lambda_n()const {
  return m_B_crf_lambda_n;
}

/**
 * Perturbation:
 * lambda value for lj interaction for state A to the power nlam-1
 */
inline double const & interaction::Perturbed_Nonbonded_Term::A_lj_lambda_n_1()const {
  return m_A_lj_lambda_n_1;
}

/**
 * Perturbation:
 * lambda value for crf interaction for state A to the power nlam-1
 */
inline double const & interaction::Perturbed_Nonbonded_Term::A_crf_lambda_n_1()const {
  return m_A_crf_lambda_n_1;
}

/**
 * Perturbation:
 * lambda value for lj interaction for state B to the power nlam-1
 */
inline double const & interaction::Perturbed_Nonbonded_Term::B_lj_lambda_n_1()const {
  return m_B_lj_lambda_n_1;
}

/**
 * Perturbation:
 * lambda value for crf interaction for state B to the power nlam-1
 */
inline double const & interaction::Perturbed_Nonbonded_Term::B_crf_lambda_n_1()const {
  return m_B_crf_lambda_n_1;
}

/**
 * Perturbation:
 * set the lambdas
 */
inline void interaction::Perturbed_Nonbonded_Term::
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
 * Perturbation:
 * lambda exponent
 */
inline int const & interaction::Perturbed_Nonbonded_Term::n()const {
  return m_n;
}

/**
 * Perturbation:
 * set the lambdas for 2d perturbations
 */
inline void interaction::Perturbed_Nonbonded_Term::
set_lambda_2d( double const lami_lj, double const lamj_lj,
            double const lamI_ljs, double const lamJ_ljs,
            double const lami_crf, double const lamj_crf,
            double const lamI_crfs, double const lamJ_crfs,
            double const lamI_lj_deriv, double const lamJ_lj_deriv,
            double const lamI_ljs_deriv, double const lamJ_ljs_deriv,
            double const lamI_crf_deriv, double const lamJ_crf_deriv,
            double const lamI_crfs_deriv, double const lamJ_crfs_deriv) {

  lamI_lj = lami_lj;
  lamJ_lj = lamj_lj;
  lamI_crf = lami_crf;
  lamJ_crf = lamj_crf;

//  std::cout << "#set_lambda_2d, lamI_lj " << std::endl;
//  std::cout << lamI_lj << " " << lamI_lj_deriv << " " << lamI_ljs_deriv << std::endl;
//  std::cout << lamJ_lj << " " << lamJ_lj_deriv << " " << lamJ_ljs_deriv << std::endl;

  // LJ
  lamA_lj = (1-lamI_lj) * (1-lamJ_lj);
  lamAB_lj = (1-lamI_lj) * lamJ_lj;
  lamBA_lj = lamI_lj * (1-lamJ_lj);
  lamB_lj = lamI_lj * lamJ_lj;

  lamA_ljs = (lamI_ljs + lamJ_ljs) * (lamI_ljs + lamJ_ljs);
  lamAB_ljs = (lamI_ljs + (1-lamJ_ljs)) * (lamI_ljs + (1-lamJ_ljs));
  lamBA_ljs = ((1-lamI_ljs) + lamJ_ljs) * ((1-lamI_ljs) + lamJ_ljs);
  lamB_ljs = ((1-lamI_ljs) + (1-lamJ_ljs)) * ((1-lamI_ljs) + (1-lamJ_ljs));

  lamA_lj_deriv = -(1-lamJ_lj)*lamI_lj_deriv - (1-lamI_lj)*lamJ_lj_deriv;
  lamAB_lj_deriv = -lamJ_lj*lamI_lj_deriv + (1-lamI_lj)*lamJ_lj_deriv;
  lamBA_lj_deriv = (1-lamJ_lj)*lamI_lj_deriv - lamI_lj* lamJ_lj_deriv;
  lamB_lj_deriv = lamJ_lj*lamI_lj_deriv + lamI_lj*lamJ_lj_deriv;

  dfA_dlamI_ljs = 2*(lamI_ljs+lamJ_ljs);
  dfAB_dlamI_ljs = 2*(1+lamI_ljs-lamJ_ljs);
  dfBA_dlamI_ljs = -2*(1-lamI_ljs+lamJ_ljs);
  dfB_dlamI_ljs = -2*(2-lamI_ljs-lamJ_ljs);

  dfA_dlamJ_ljs = dfA_dlamI_ljs;
  dfAB_dlamJ_ljs = -dfAB_dlamI_ljs;
  dfBA_dlamJ_ljs = -dfBA_dlamI_ljs;
  dfB_dlamJ_ljs = dfB_dlamI_ljs;

  lamA_ljs_deriv = dfA_dlamI_ljs * lamI_ljs_deriv + dfA_dlamJ_ljs * lamJ_ljs_deriv;
  lamAB_ljs_deriv = dfAB_dlamI_ljs * lamI_ljs_deriv + dfAB_dlamJ_ljs * lamJ_ljs_deriv;
  lamBA_ljs_deriv = dfBA_dlamI_ljs * lamI_ljs_deriv + dfBA_dlamJ_ljs * lamJ_ljs_deriv;
  lamB_ljs_deriv = dfB_dlamI_ljs * lamI_ljs_deriv + dfB_dlamJ_ljs * lamJ_ljs_deriv;

  // CRF
  lamA_crf = (1-lamI_crf) * (1-lamJ_crf);
  lamAB_crf = (1-lamI_crf) * lamJ_crf;
  lamBA_crf = lamI_crf * (1-lamJ_crf);
  lamB_crf = lamI_crf * lamJ_crf;

  lamA_crfs = (lamI_crfs + lamJ_crfs) * (lamI_crfs + lamJ_crfs);
  lamAB_crfs = (lamI_crfs + (1-lamJ_crfs)) * (lamI_crfs + (1-lamJ_crfs));
  lamBA_crfs = ((1-lamI_crfs) + lamJ_crfs) * ((1-lamI_crfs) + lamJ_crfs);
  lamB_crfs = ((1-lamI_crfs) + (1-lamJ_crfs)) * ((1-lamI_crfs) + (1-lamJ_crfs));

  lamA_crf_deriv = -(1-lamJ_crf)*lamI_crf_deriv - (1-lamI_crf)*lamJ_crf_deriv;
  lamAB_crf_deriv = -lamJ_crf*lamI_crf_deriv + (1-lamI_crf)*lamJ_crf_deriv;
  lamBA_crf_deriv = (1-lamJ_crf)*lamI_crf_deriv - lamI_crf* lamJ_crf_deriv;
  lamB_crf_deriv = lamJ_crf*lamI_crf_deriv + lamI_crf*lamJ_crf_deriv;

  dfA_dlamI_crfs = 2*(lamI_crfs+lamJ_crfs);
  dfAB_dlamI_crfs = 2*(1+lamI_crfs-lamJ_crfs);
  dfBA_dlamI_crfs = -2*(1-lamI_crfs+lamJ_crfs);
  dfB_dlamI_crfs = -2*(2-lamI_crfs-lamJ_crfs);

  dfA_dlamJ_crfs = dfA_dlamI_crfs;
  dfAB_dlamJ_crfs = -dfAB_dlamI_crfs;
  dfBA_dlamJ_crfs = -dfBA_dlamI_crfs;
  dfB_dlamJ_crfs = dfB_dlamI_crfs;

  lamA_crfs_deriv = dfA_dlamI_crfs * lamI_crfs_deriv + dfA_dlamJ_crfs * lamJ_crfs_deriv;
  lamAB_crfs_deriv = dfAB_dlamI_crfs * lamI_crfs_deriv + dfAB_dlamJ_crfs * lamJ_crfs_deriv;
  lamBA_crfs_deriv = dfBA_dlamI_crfs * lamI_crfs_deriv + dfBA_dlamJ_crfs * lamJ_crfs_deriv;
  lamB_crfs_deriv = dfB_dlamI_crfs * lamI_crfs_deriv + dfB_dlamJ_crfs * lamJ_crfs_deriv;

}

/**
 * Perturbation:
 * set the lambdas
 */
inline void interaction::Perturbed_Nonbonded_Term::
set_lambda_2d_other(double const LJ_lambda, double const LJS_lambda,
        double const CRF_lambda, double const CRFS_lambda,
        double const lj_lambda_derivative, double const ljs_lambda_derivative,
        double const crf_lambda_derivative, double const crfs_lambda_derivative,
        int const n) {

  lj_lambda = LJ_lambda;
  ljs_lambda = LJS_lambda;
  crf_lambda = CRF_lambda;
  crfs_lambda = CRFS_lambda;

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
 * the reaction field contribution for a given pair
 * using the soft interaction
 */
inline void interaction::Perturbed_Nonbonded_Term
::rf_soft_interaction(math::Vec const &r, double const A_q, double const B_q,
        double const alpha_crf, math::Vec & force, double &e_rf, double & de_rf,
        bool selfterm_correction, unsigned int eps) {
  const double dist2 = abs2(r);

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

  force = (m_A_crf_lambda_n * A_q * A_crf_cut3i +
          m_B_crf_lambda_n * B_q * B_crf_cut3i) * math::four_pi_eps_i *r;

  double const A_e_rf = A_q * (-A_crf_2cut3i * dist2 - m_crf_cut[eps]);
  double const B_e_rf = B_q * (-B_crf_2cut3i * dist2 - m_crf_cut[eps]);

  e_rf = (m_A_crf_lambda_n * A_e_rf + m_B_crf_lambda_n * B_e_rf) * math::four_pi_eps_i;

  if (selfterm_correction)
    e_rf += A_q * math::four_pi_eps_i * m_crf_cut[eps];
  // Chris: CHECK! I'm not sure if the self-term correction is not wrong... 
  //        There is no B state info going in here at all!
  //        If a B state is added, then there should also be a contribution to de_rf

  DEBUG(11, "A_crf_pert: " << A_crf_pert << " B_crf_pert: " << B_crf_pert);
  DEBUG(11, "m_lambda_exp: " << m_lambda_exp);

  de_rf = ((m_A_crf_lambda_n * A_q * m_B_crfs_lambda * A_crf_pert -
          m_B_crf_lambda_n * B_q * m_A_crfs_lambda * B_crf_pert) * dist2 * alpha_crf +
          (m_B_crf_lambda_n_1 * B_e_rf -
          m_A_crf_lambda_n_1 * A_e_rf) * m_lambda_exp) * math::four_pi_eps_i;
}
/**
 * helper function to calculate the force and energy for 
 * the reaction field contribution for a given pair
 * using the soft interaction (for extendedTI) 
 */
inline void interaction::Perturbed_Nonbonded_Term
::rf_soft_interaction_ext(math::Vec const &r, double const A_q, double const B_q,
        double const alpha_crf, double & A_e_rf, double & B_e_rf,
        double & A_de_rf, double & B_de_rf, double const crfs_lambda,
        unsigned int eps) {
  const double dist2 = abs2(r);

  const double A_cut2soft = m_cut2 + alpha_crf * m_B_crfs_lambda2;
  const double B_cut2soft = m_cut2 + alpha_crf * m_A_crfs_lambda2;

  const double A_cut2soft3 = A_cut2soft * A_cut2soft * A_cut2soft;
  const double B_cut2soft3 = B_cut2soft * B_cut2soft * B_cut2soft;

  const double A_crf_2cut3i = m_crf_2[eps] / sqrt(A_cut2soft3);
  const double B_crf_2cut3i = m_crf_2[eps] / sqrt(B_cut2soft3);

  const double A_crf_pert = 3.0 * A_crf_2cut3i / A_cut2soft;
  const double B_crf_pert = 3.0 * B_crf_2cut3i / B_cut2soft;

  A_e_rf = A_q * (-A_crf_2cut3i * dist2 - m_crf_cut[eps]) * math::four_pi_eps_i;
  B_e_rf = B_q * (-B_crf_2cut3i * dist2 - m_crf_cut[eps]) * math::four_pi_eps_i;

  A_de_rf = A_q * crfs_lambda * A_crf_pert * dist2 *
             math::four_pi_eps_i * alpha_crf;

  B_de_rf = B_q * (1 - crfs_lambda) * B_crf_pert * dist2 *
             math::four_pi_eps_i * alpha_crf;

}

/**
 * helper function to calculate the force and energy for 
 * the reaction field contribution for a given pair
 * using the soft interaction, for 2d perturbations
 */
inline void interaction::Perturbed_Nonbonded_Term
::rf_soft_2d_interaction(math::Vec const &r, double const A_q, double const B_q,
        double const alpha_crf, math::Vec & force, double &e_rf, double & de_rf,
        double &pert2d_crf, double &pert2d_crfs,
        bool selfterm_correction, unsigned int eps) {
  const double dist2 = abs2(r);

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

  force = (m_A_crf_lambda_n * A_q * A_crf_cut3i +
          m_B_crf_lambda_n * B_q * B_crf_cut3i) * math::four_pi_eps_i *r;

  double const A_e_rf = A_q * (-A_crf_2cut3i * dist2 - m_crf_cut[eps]);
  double const B_e_rf = B_q * (-B_crf_2cut3i * dist2 - m_crf_cut[eps]);

  e_rf = (m_A_crf_lambda_n * A_e_rf + m_B_crf_lambda_n * B_e_rf) * math::four_pi_eps_i;

  if (selfterm_correction)
    e_rf += A_q * math::four_pi_eps_i * m_crf_cut[eps];
  // Chris: CHECK! I'm not sure if the self-term correction is not wrong... 
  //        There is no B state info going in here at all!
  //        If a B state is added, then there should also be a contribution to de_rf

  DEBUG(11, "A_crf_pert: " << A_crf_pert << " B_crf_pert: " << B_crf_pert);
  DEBUG(11, "m_lambda_exp: " << m_lambda_exp);

  de_rf = ((m_A_crf_lambda_n * A_q * m_B_crfs_lambda * A_crf_pert -
          m_B_crf_lambda_n * B_q * m_A_crfs_lambda * B_crf_pert) * dist2 * alpha_crf +
          (m_B_crf_lambda_n_1 * B_e_rf -
          m_A_crf_lambda_n_1 * A_e_rf) * m_lambda_exp) * math::four_pi_eps_i;

  pert2d_crf = (pow(crf_lambda, m_lambda_exp-1) * B_e_rf - pow(1-crf_lambda, m_lambda_exp-1) * A_e_rf)
               * m_lambda_exp * math::four_pi_eps_i;
  pert2d_crfs = (m_A_crf_lambda_n * A_q * crfs_lambda * A_crf_pert -
                m_B_crf_lambda_n * B_q * (1-crfs_lambda) * B_crf_pert) * dist2 * alpha_crf *
                math::four_pi_eps_i;
}

inline void interaction::Perturbed_Nonbonded_Term      
::pol_rf_self_soft_interaction(double const A_qi,double const B_qi,
        double & e_rf, double & de_rf, bool selfterm_correction ,unsigned int eps){ 

  double const field_term = math::four_pi_eps_i * (-m_crf_cut[eps]);
  double const charge_term = (m_A_crfs_lambda * A_qi) + (m_B_crfs_lambda * B_qi);
 

  e_rf =  math::four_pi_eps_i * charge_term * charge_term * (-m_crf_cut[eps]); 

  if (selfterm_correction) 
    e_rf += A_qi*A_qi*math::four_pi_eps_i * m_crf_cut[eps];          
  
  de_rf = field_term * charge_term * (B_qi - A_qi) *2;      

}

/**
 * helper function to calculate the force and energy for 
 * the reaction field contribution for a given pair
 * using the soft interaction (with polarisation)
 */
inline void interaction::Perturbed_Nonbonded_Term
::pol_rf_soft_interaction(math::Vec const &r,
        math::Vec const &rp1,
        math::Vec const &rp2,
        math::Vec const &rpp,
        double const A_qi, double const A_qj,
        double const B_qi, double const B_qj,
        double cqi, double cqj,
        double const alpha_crf,
        double force[], double &e_rf, double & de_rf,
        bool selfterm_correction, unsigned int eps) {
  const double dist2 = abs2(r);
  const double dist2p1 = abs2(rp1);
  const double dist2p2 = abs2(rp2);
  const double dist2pp = abs2(rpp);

  const double A_qeps = (A_qi - cqi) * (A_qj - cqj) * math::four_pi_eps_i;
  const double A_qepsp1 = (A_qi - cqi) * cqj * math::four_pi_eps_i;
  const double A_qepsp2 = (A_qj - cqj) * cqi * math::four_pi_eps_i;
  const double A_qepspp = cqi * cqj * math::four_pi_eps_i;

  const double B_qeps = (B_qi - cqi) * (B_qj - cqj) * math::four_pi_eps_i;
  const double B_qepsp1 = (B_qi - cqi) * cqj * math::four_pi_eps_i;
  const double B_qepsp2 = (B_qj - cqj) * cqi * math::four_pi_eps_i;
  const double B_qepspp = cqi * cqj * math::four_pi_eps_i;

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

  const double A_lambda_cut = m_A_crf_lambda_n * A_crf_cut3i;
  const double B_lambda_cut = m_B_crf_lambda_n * B_crf_cut3i;

  force[0] = A_qeps * A_lambda_cut + B_qeps*B_lambda_cut;
  force[1] = A_qepsp1 * A_lambda_cut + B_qepsp1*B_lambda_cut;
  force[2] = A_qepsp2 * A_lambda_cut + B_qepsp2*B_lambda_cut;
  force[3] = A_qepspp * A_lambda_cut + B_qepspp*B_lambda_cut;

  double const A_erf = A_qeps * (-A_crf_2cut3i * dist2 - m_crf_cut[eps]);
  double const A_erfp1 = A_qepsp1 * (-A_crf_2cut3i * dist2p1 - m_crf_cut[eps]);
  double const A_erfp2 = A_qepsp2 * (-A_crf_2cut3i * dist2p2 - m_crf_cut[eps]);
  double const A_erfpp = A_qepspp * (-A_crf_2cut3i * dist2pp - m_crf_cut[eps]);

  double const B_erf = B_qeps * (-B_crf_2cut3i * dist2 - m_crf_cut[eps]);
  double const B_erfp1 = B_qepsp1 * (-B_crf_2cut3i * dist2p1 - m_crf_cut[eps]);
  double const B_erfp2 = B_qepsp2 * (-B_crf_2cut3i * dist2p2 - m_crf_cut[eps]);
  double const B_erfpp = B_qepspp * (-B_crf_2cut3i * dist2pp - m_crf_cut[eps]);

  e_rf = m_A_crf_lambda_n * (A_erf + A_erfp1 + A_erfp2 + A_erfpp)
          + m_B_crf_lambda_n * (B_erf + B_erfp1 + B_erfp2 + B_erfpp);

  if (selfterm_correction)
    e_rf += A_qeps * m_crf_cut[eps];

  de_rf = (m_A_crf_lambda_n * A_qeps * m_B_crfs_lambda * A_crf_pert -
          m_B_crf_lambda_n * B_qeps * m_A_crfs_lambda * B_crf_pert) * dist2 * alpha_crf +
          (m_B_crf_lambda_n_1 * B_erf - m_A_crf_lambda_n_1 * A_erf) * m_lambda_exp
          + (m_A_crf_lambda_n * A_qepsp1 * m_B_crfs_lambda * A_crf_pert -
          m_B_crf_lambda_n * B_qepsp1 * m_A_crfs_lambda * B_crf_pert) * dist2p1 * alpha_crf +
          (m_B_crf_lambda_n_1 * B_erfp1 - m_A_crf_lambda_n_1 * A_erfp1) * m_lambda_exp
          + (m_A_crf_lambda_n * A_qepsp2 * m_B_crfs_lambda * A_crf_pert -
          m_B_crf_lambda_n * B_qepsp2 * m_A_crfs_lambda * B_crf_pert) * dist2p2 * alpha_crf +
          (m_B_crf_lambda_n_1 * B_erfp2 - m_A_crf_lambda_n_1 * A_erfp2) * m_lambda_exp
          + (m_A_crf_lambda_n * A_qepspp * m_B_crfs_lambda * A_crf_pert -
          m_B_crf_lambda_n * B_qepspp * m_A_crfs_lambda * B_crf_pert) * dist2pp * alpha_crf +
          (m_B_crf_lambda_n_1 * B_erfpp - m_A_crf_lambda_n_1 * A_erfpp) * m_lambda_exp;

}

inline void interaction::Perturbed_Nonbonded_Term
::rf_scaled_interaction(math::Vec const &r, double const A_q, double const B_q,
        double const l, double const alpha_crf,
        double const A_scale, double const B_scale,
        math::Vec & force, double &e_rf, double & de_rf,
        bool selfterm_correction, unsigned int eps) {
  const double dist2 = abs2(r);

  const double A_cut2soft = m_cut2 + alpha_crf * m_B_crfs_lambda2;
  const double B_cut2soft = m_cut2 + alpha_crf * m_A_crfs_lambda2;

  const double A_cut2soft3 = A_cut2soft * A_cut2soft*A_cut2soft;
  const double B_cut2soft3 = B_cut2soft * B_cut2soft*B_cut2soft;

  const double A_crf_2cut3i = m_crf_2[eps] / sqrt(A_cut2soft3);
  const double B_crf_2cut3i = m_crf_2[eps] / sqrt(B_cut2soft3);

  const double A_crf_cut3i = 2 * A_crf_2cut3i;
  const double B_crf_cut3i = 2 * B_crf_2cut3i;

  const double A_crf_pert = 3.0 * A_crf_2cut3i / A_cut2soft;
  const double B_crf_pert = 3.0 * B_crf_2cut3i / B_cut2soft;

  const double A_scaled_lambda_n = A_scale * m_A_crf_lambda_n;
  const double B_scaled_lambda_n = B_scale * m_B_crf_lambda_n;

  force = (A_scaled_lambda_n * A_q * A_crf_cut3i +
          B_scaled_lambda_n * B_q * B_crf_cut3i) * math::four_pi_eps_i *r;

  double const A_e_rf = A_q * (-A_crf_2cut3i * dist2 - m_crf_cut[eps]);
  double const B_e_rf = B_q * (-B_crf_2cut3i * dist2 - m_crf_cut[eps]);

  e_rf = (A_scaled_lambda_n * A_e_rf + B_scaled_lambda_n * B_e_rf) * math::four_pi_eps_i;
  if (selfterm_correction)
    e_rf += A_q * math::four_pi_eps_i * m_crf_cut[eps];

  de_rf = ((A_scaled_lambda_n * A_q * m_B_crfs_lambda * A_crf_pert -
          B_scaled_lambda_n * B_q * m_A_crfs_lambda * B_crf_pert) * dist2 * alpha_crf +
          (B_scale * m_B_crf_lambda_n_1 * B_e_rf -
          A_scale * m_A_crf_lambda_n_1 * A_e_rf) * m_lambda_exp) * math::four_pi_eps_i;

}

inline void interaction::Perturbed_Nonbonded_Term
::lj_crf_soft_interaction(math::Vec const &r,
        double const A_c6, double const A_c12,
        double const B_c6, double const B_c12,
        double const A_q, double const B_q,
        double const alpha_lj, double const alpha_crf,
        double &force1, double &force6, double &force12,
        double &e_lj, double &e_crf,
        double &de_lj, double & de_crf, unsigned int eps) {
  double A_c126, B_c126;

//  std::cout << "\tA_c6 " << A_c6 << " A_c12 " << A_c12 << " B_c6 " << B_c6 << " B_c12 " << B_c12 << std::endl;

  if (A_c6 != 0) A_c126 = A_c12 / A_c6;
  else A_c126 = 0.0;
  if (B_c6 != 0) B_c126 = B_c12 / B_c6;
  else B_c126 = 0.0;

  const double dist2 = abs2(r);
  assert(dist2 != 0);

  const double A_dist2soft = dist2 + alpha_crf*m_B_crfs_lambda2;
  const double B_dist2soft = dist2 + alpha_crf*m_A_crfs_lambda2;

  const double A_distisoft = 1.0 / sqrt(A_dist2soft);
  const double B_distisoft = 1.0 / sqrt(B_dist2soft);

  const double A_dist3isoft = A_distisoft / A_dist2soft;
  const double B_dist3isoft = B_distisoft / B_dist2soft;

  const double dist4 = dist2 * dist2;
  const double dist6 = dist4 * dist2;

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

  // substitute A_dist3isoft thing. just like here -- daniel
  force1 = (m_A_crf_lambda_n * A_q * (A_dist3isoft + A_crf_cut3i) +
          m_B_crf_lambda_n * B_q * (B_dist3isoft + B_crf_cut3i)) *
          math::four_pi_eps_i;

  force6 = -6.0 * (m_A_lj_lambda_n * A_c6 * A_dist6isoft * A_dist6isoft +
          m_B_lj_lambda_n * B_c6 * B_dist6isoft * B_dist6isoft) * dist4;

  force12 = 12 * (m_A_lj_lambda_n * A_c12 * A_dist6isoft * A_dist6isoft * A_dist6isoft +
          m_B_lj_lambda_n * B_c12 * B_dist6isoft * B_dist6isoft * B_dist6isoft) * dist4;

  const double A_e_lj = (A_c12 * A_dist6isoft - A_c6) * A_dist6isoft;
  const double B_e_lj = (B_c12 * B_dist6isoft - B_c6) * B_dist6isoft;

  const double A_e_crf = A_q * (A_distisoft - A_crf_2cut3i * dist2 - m_crf_cut[eps]);
  const double B_e_crf = B_q * (B_distisoft - B_crf_2cut3i * dist2 - m_crf_cut[eps]);

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

//  std::cout << "\tA_e_lj " <<  A_e_lj << " m_A_lj_lambda_n " <<  m_A_lj_lambda_n <<
//               " B_e_lj " << B_e_lj << " m_B_lj_lambda_n " << m_B_lj_lambda_n << std::endl;

  e_lj = m_A_lj_lambda_n * A_e_lj + m_B_lj_lambda_n * B_e_lj;

  e_crf = (m_A_crf_lambda_n * A_e_crf + m_B_crf_lambda_n * B_e_crf) * math::four_pi_eps_i;

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

//  pert2d_lj = m_lambda_exp * (pow(lj_lambda,n-1) * B_e_lj - pow(1-lj_lambda,n-1) * A_e_lj);
//  pert2d_ljs = -2.0 * alpha_lj * (m_A_lj_lambda_n * ljs_lambda * A_c126 * A_dist6isoft * A_dist6isoft *
//          (2 * A_c12 * A_dist6isoft - A_c6) -
//          m_B_lj_lambda_n * (1-ljs_lambda) * B_c126 * B_dist6isoft * B_dist6isoft *
//          (2 * B_c12 * B_dist6isoft - B_c6));
   DEBUG(15, "nonbonded energy = " << (e_lj+e_crf) << ", force = " << (force1+force6+force12));

}

/* ############## */
/* for extendedTI */
/* ############## */

inline void interaction::Perturbed_Nonbonded_Term
::lj_crf_soft_interaction_ext(math::Vec const &r,
        double const A_c6, double const A_c12,
        double const B_c6, double const B_c12,
        double const A_q, double const B_q,
        double const alpha_lj, double const alpha_crf,
        double &A_e_lj, double &B_e_lj,
        double &A_e_crf, double &B_e_crf,
        double &A_de_lj, double &B_de_lj,
        double &A_de_crf, double &B_de_crf,
        double const lam, unsigned int eps) {

  const double ljs_lambda = lam;
  const double A_ljs_lambda2 = (1 - ljs_lambda) * (1 - ljs_lambda);
  const double B_ljs_lambda2 = ljs_lambda * ljs_lambda;

  const double crfs_lambda = lam;
  const double A_crfs_lambda2 = (1 - crfs_lambda) * (1 - crfs_lambda);
  const double B_crfs_lambda2 = crfs_lambda * crfs_lambda;

  double A_c126, B_c126;

  if (A_c6 != 0) A_c126 = A_c12 / A_c6;
  else A_c126 = 0.0;
  if (B_c6 != 0) B_c126 = B_c12 / B_c6;
  else B_c126 = 0.0;

  const double dist2 = abs2(r);
  assert(dist2 != 0);

  const double A_dist2soft = dist2 + alpha_crf*B_crfs_lambda2;
  const double B_dist2soft = dist2 + alpha_crf*A_crfs_lambda2;

  const double A_distisoft = 1.0 / sqrt(A_dist2soft);
  const double B_distisoft = 1.0 / sqrt(B_dist2soft);

  const double A_dist3isoft = A_distisoft / A_dist2soft;
  const double B_dist3isoft = B_distisoft / B_dist2soft;

  const double dist4 = dist2 * dist2;
  const double dist6 = dist4 * dist2;

  const double A_dist6soft = dist6 + alpha_lj * B_ljs_lambda2*A_c126;
  const double B_dist6soft = dist6 + alpha_lj * A_ljs_lambda2*B_c126;

  const double A_dist6isoft = 1.0 / A_dist6soft;
  const double B_dist6isoft = 1.0 / B_dist6soft;

  const double A_cut2soft = m_cut2 + alpha_crf * B_crfs_lambda2;
  const double B_cut2soft = m_cut2 + alpha_crf * A_crfs_lambda2;

  const double A_cut2soft3 = A_cut2soft * A_cut2soft * A_cut2soft;
  const double B_cut2soft3 = B_cut2soft * B_cut2soft * B_cut2soft;

  const double A_crf_2cut3i = m_crf_2[eps] / sqrt(A_cut2soft3);
  const double B_crf_2cut3i = m_crf_2[eps] / sqrt(B_cut2soft3);


  const double A_crf_pert = 3.0 * A_crf_2cut3i / A_cut2soft;
  const double B_crf_pert = 3.0 * B_crf_2cut3i / B_cut2soft;

  DEBUG(10, "just checking some values \n  A_ljs_lambda2 " << A_ljs_lambda2
          << "\n  B_ljs_lambda2 " << B_ljs_lambda2
          << "\n  A_crfs_lambda2 " << A_crfs_lambda2
          << "\n  B_crfs_lambda2 " << B_crfs_lambda2
          << "\n  A_dist6isoft " << A_dist6isoft
          << "\n  B_dist6isoft " << B_dist6isoft);

  A_e_lj = (A_c12 * A_dist6isoft - A_c6) * A_dist6isoft;
  B_e_lj = (B_c12 * B_dist6isoft - B_c6) * B_dist6isoft;

  A_e_crf = A_q * (A_distisoft - A_crf_2cut3i * dist2 - m_crf_cut[eps]) * math::four_pi_eps_i;
  B_e_crf = B_q * (B_distisoft - B_crf_2cut3i * dist2 - m_crf_cut[eps]) * math::four_pi_eps_i;

  A_de_lj = -2.0 * alpha_lj * ljs_lambda * A_c126 * A_dist6isoft * A_dist6isoft *
            (2 * A_c12 * A_dist6isoft - A_c6);

  B_de_lj = -2.0 * alpha_lj * (1 - ljs_lambda) * B_c126 * B_dist6isoft * B_dist6isoft *
            (2 * B_c12 * B_dist6isoft - B_c6);

  A_de_crf = -A_q * crfs_lambda * (A_dist3isoft - A_crf_pert * dist2) *
             math::four_pi_eps_i * alpha_crf;

  B_de_crf = -B_q * (1 - crfs_lambda) * (B_dist3isoft - B_crf_pert * dist2) *
             math::four_pi_eps_i * alpha_crf;

  DEBUG(10, "just checking\nA_e_lj " << A_e_lj
          << "\nB_e_lj " << B_e_lj
          << "\nA_e_crf " << A_e_crf
          << "\nB_e_crf " << B_e_crf
          << "\nA_de_lj " << A_de_lj
          << "\nB_de_lj " << B_de_lj
          << "\nA_de_crf " << A_de_crf
          << "\nB_de_crf " << B_de_crf);
  
} // extendedTI

// 2D perturbation
inline void interaction::Perturbed_Nonbonded_Term
::lj_crf_soft_2d_interaction(math::Vec const &r,
        double const A_i_c6, double const A_i_c12,
        double const B_i_c6, double const B_i_c12,
        double const A_j_c6, double const A_j_c12,
        double const B_j_c6, double const B_j_c12,
        double const A_i_q, double const B_i_q,
        double const A_j_q, double const B_j_q,
        double const alpha_lj, double const alpha_crf,
        double &force1, double &force6, double &force12,
        double &e_lj, double &e_crf,
        double &de_lj, double & de_crf, 
        double &pert2d_lamI_lj, double &pert2d_lamI_ljs, 
        double &pert2d_lamI_crf, double &pert2d_lamI_crfs,
        double &pert2d_lamJ_lj, double &pert2d_lamJ_ljs,
        double &pert2d_lamJ_crf, double &pert2d_lamJ_crfs, unsigned int eps) {

  const double dist2 = abs2(r);
  assert(dist2 != 0);

  //LJ
  const double A_c6 = sqrt(A_i_c6 * A_j_c6);
  const double A_c12 = sqrt(A_i_c12 * A_j_c12);

  const double AB_c6 = sqrt(A_i_c6 * B_j_c6);
  const double AB_c12 = sqrt(A_i_c12 * B_j_c12);

  const double BA_c6 = sqrt(B_i_c6 * A_j_c6);
  const double BA_c12 = sqrt(B_i_c12 * A_j_c12);

  const double B_c6 = sqrt(B_i_c6 * B_j_c6);
  const double B_c12 = sqrt(B_i_c12 * B_j_c12);

//  std::cout << "\tA_i_c6 " << A_i_c6 << " A_j_c6 " << A_j_c6 << std::endl;
//  std::cout << "\tA_c6 " << A_c6 << " A_c12 " << A_c12 << " B_c6 " << B_c6 << " B_c12 " << B_c12 << std::endl;

  double A_c126, AB_c126, BA_c126, B_c126;
  if (A_c6 != 0) A_c126 = A_c12 / A_c6;
  else A_c126 = 0.0;
  if (AB_c6 != 0) AB_c126 = AB_c12 / AB_c6;
  else AB_c126 = 0.0;
  if (BA_c6 != 0) BA_c126 = BA_c12 / BA_c6;
  else BA_c126 = 0.0;
  if (B_c6 != 0) B_c126 = B_c12 / B_c6;
  else B_c126 = 0.0;

  const double dist4 = dist2 * dist2;
  const double dist6 = dist4 * dist2;

  const double A_dist6soft = dist6 + alpha_lj * lamA_ljs*A_c126;
  const double AB_dist6soft= dist6 + alpha_lj * lamAB_ljs*AB_c126;
  const double BA_dist6soft= dist6 + alpha_lj * lamBA_ljs*BA_c126;
  const double B_dist6soft = dist6 + alpha_lj * lamB_ljs*B_c126;

  const double A_dist6isoft = 1.0 / A_dist6soft;
  const double AB_dist6isoft = 1.0 / AB_dist6soft;
  const double BA_dist6isoft = 1.0 / BA_dist6soft;
  const double B_dist6isoft = 1.0 / B_dist6soft;

  force6 = -6.0 * (lamA_lj * A_c6 * A_dist6isoft * A_dist6isoft +
          lamAB_lj * AB_c6 * AB_dist6isoft * AB_dist6isoft +
          lamBA_lj * BA_c6 * BA_dist6isoft * BA_dist6isoft +
          lamB_lj * B_c6 * B_dist6isoft * B_dist6isoft) * dist4;

  force12 = 12 * (lamA_lj * A_c12 * A_dist6isoft * A_dist6isoft * A_dist6isoft +
          lamAB_lj * AB_c12 * AB_dist6isoft * AB_dist6isoft * AB_dist6isoft +
          lamBA_lj * BA_c12 * BA_dist6isoft * BA_dist6isoft * BA_dist6isoft +
          lamB_lj * B_c12 * B_dist6isoft * B_dist6isoft * B_dist6isoft) * dist4;

  const double A_e_lj = (A_c12 * A_dist6isoft - A_c6) * A_dist6isoft;
  const double B_e_lj = (B_c12 * B_dist6isoft - B_c6) * B_dist6isoft;
  const double AB_e_lj = (AB_c12 * AB_dist6isoft - AB_c6) * AB_dist6isoft;
  const double BA_e_lj = (BA_c12 * BA_dist6isoft - BA_c6) * BA_dist6isoft;

//  std::cout << "\tA_e_lj " <<  A_e_lj << "  lamA_lj " <<  lamA_lj << 
//               " AB_e_lj " << AB_e_lj << " lamAB_lj " << lamAB_lj <<
//               " BA_e_lj " << BA_e_lj << " lamBA_lj " << lamBA_lj <<
//               "  B_e_lj " <<  B_e_lj << "  lamB_lj " <<  lamB_lj << std::endl;
  e_lj = lamA_lj * A_e_lj + lamAB_lj * AB_e_lj +
         lamBA_lj * BA_e_lj + lamB_lj * B_e_lj;

  de_lj = A_e_lj * lamA_lj_deriv - (alpha_lj * lamA_lj * A_c126 * A_dist6isoft * A_dist6isoft) *
          (2 * A_c12 * A_dist6isoft - A_c6) * lamA_ljs_deriv +
          AB_e_lj * lamAB_lj_deriv - (alpha_lj * lamAB_lj * AB_c126 * AB_dist6isoft * AB_dist6isoft) *
          (2 * AB_c12 * AB_dist6isoft - AB_c6) * lamAB_ljs_deriv +
          BA_e_lj * lamBA_lj_deriv - (alpha_lj * lamBA_lj * BA_c126 * BA_dist6isoft * BA_dist6isoft) *
          (2 * BA_c12 * BA_dist6isoft - BA_c6) * lamBA_ljs_deriv +
          B_e_lj * lamB_lj_deriv - (alpha_lj * lamB_lj * B_c126 * B_dist6isoft * B_dist6isoft) *
          (2 * B_c12 * B_dist6isoft - B_c6) * lamB_ljs_deriv;

  pert2d_lamI_lj  = A_e_lj * -(1 - lamJ_lj) + AB_e_lj * -lamJ_lj +
                 BA_e_lj * (1 - lamJ_lj) + B_e_lj * lamJ_lj; 

  pert2d_lamI_ljs = - (alpha_lj * lamA_lj * A_c126 * A_dist6isoft * A_dist6isoft) *
          (2 * A_c12 * A_dist6isoft - A_c6) * dfA_dlamI_ljs 
           - (alpha_lj * lamAB_lj * AB_c126 * AB_dist6isoft * AB_dist6isoft) *
          (2 * AB_c12 * AB_dist6isoft - AB_c6) * dfAB_dlamI_ljs  
          - (alpha_lj * lamBA_lj * BA_c126 * BA_dist6isoft * BA_dist6isoft) *
          (2 * BA_c12 * BA_dist6isoft - BA_c6) * dfBA_dlamI_ljs 
          - (alpha_lj * lamB_lj * B_c126 * B_dist6isoft * B_dist6isoft) *
          (2 * B_c12 * B_dist6isoft - B_c6) * dfB_dlamI_ljs;

  pert2d_lamJ_lj  = A_e_lj * -(1 - lamI_lj) + AB_e_lj * (1 - lamI_lj) +
                 BA_e_lj * -lamI_lj + B_e_lj * lamI_lj; 

  pert2d_lamJ_ljs = - (alpha_lj * lamA_lj * A_c126 * A_dist6isoft * A_dist6isoft) *
          (2 * A_c12 * A_dist6isoft - A_c6) * dfA_dlamJ_ljs 
           - (alpha_lj * lamAB_lj * AB_c126 * AB_dist6isoft * AB_dist6isoft) *
          (2 * AB_c12 * AB_dist6isoft - AB_c6) * dfAB_dlamJ_ljs  
          - (alpha_lj * lamBA_lj * BA_c126 * BA_dist6isoft * BA_dist6isoft) *
          (2 * BA_c12 * BA_dist6isoft - BA_c6) * dfBA_dlamJ_ljs 
          - (alpha_lj * lamB_lj * B_c126 * B_dist6isoft * B_dist6isoft) *
          (2 * B_c12 * B_dist6isoft - B_c6) * dfB_dlamJ_ljs;

  //CRF
  double A_q = A_i_q * A_j_q;
  double AB_q = A_i_q * B_j_q;
  double BA_q = B_i_q * A_j_q;
  double B_q = B_i_q * B_j_q;

  const double A_dist2soft = dist2 + alpha_crf*lamA_crfs;
  const double AB_dist2soft = dist2 + alpha_crf*lamAB_crfs;
  const double BA_dist2soft = dist2 + alpha_crf*lamBA_crfs;
  const double B_dist2soft = dist2 + alpha_crf*lamB_crfs;

  const double A_distisoft = 1.0 / sqrt(A_dist2soft);
  const double AB_distisoft = 1.0 / sqrt(AB_dist2soft);
  const double BA_distisoft = 1.0 / sqrt(BA_dist2soft);
  const double B_distisoft = 1.0 / sqrt(B_dist2soft);

  const double A_dist3isoft = A_distisoft / A_dist2soft;
  const double AB_dist3isoft = AB_distisoft / AB_dist2soft;
  const double BA_dist3isoft = BA_distisoft / BA_dist2soft;
  const double B_dist3isoft = B_distisoft / B_dist2soft;

  const double A_cut2soft = m_cut2 + alpha_crf * lamA_crfs;
  const double AB_cut2soft = m_cut2 + alpha_crf * lamAB_crfs;
  const double BA_cut2soft = m_cut2 + alpha_crf * lamBA_crfs;
  const double B_cut2soft = m_cut2 + alpha_crf * lamB_crfs;

  const double A_cut2soft3 = A_cut2soft * A_cut2soft * A_cut2soft;
  const double AB_cut2soft3 = AB_cut2soft * AB_cut2soft * AB_cut2soft;
  const double BA_cut2soft3 = BA_cut2soft * BA_cut2soft * BA_cut2soft;
  const double B_cut2soft3 = B_cut2soft * B_cut2soft * B_cut2soft;

  const double A_crf_2cut3i = m_crf_2[eps] / sqrt(A_cut2soft3);
  const double AB_crf_2cut3i = m_crf_2[eps] / sqrt(AB_cut2soft3);
  const double BA_crf_2cut3i = m_crf_2[eps] / sqrt(BA_cut2soft3);
  const double B_crf_2cut3i = m_crf_2[eps] / sqrt(B_cut2soft3);

  const double A_crf_cut3i = 2 * A_crf_2cut3i;
  const double AB_crf_cut3i = 2 * AB_crf_2cut3i;
  const double BA_crf_cut3i = 2 * BA_crf_2cut3i;
  const double B_crf_cut3i = 2 * B_crf_2cut3i;

  const double A_crf_pert = 3.0 * A_crf_2cut3i / A_cut2soft;
  const double AB_crf_pert = 3.0 * AB_crf_2cut3i / AB_cut2soft;
  const double BA_crf_pert = 3.0 * BA_crf_2cut3i / BA_cut2soft;
  const double B_crf_pert = 3.0 * B_crf_2cut3i / B_cut2soft;

  force1 = (lamA_crf * A_q * (A_dist3isoft + A_crf_cut3i) +
          lamAB_crf * AB_q * (AB_dist3isoft + AB_crf_cut3i) +
          lamBA_crf * BA_q * (BA_dist3isoft + BA_crf_cut3i) +
          lamB_crf * B_q * (B_dist3isoft + B_crf_cut3i)) *
          math::four_pi_eps_i;

  const double A_e_crf = A_q * (A_distisoft - A_crf_2cut3i * dist2 - m_crf_cut[eps]);
  const double AB_e_crf = AB_q * (AB_distisoft - AB_crf_2cut3i * dist2 - m_crf_cut[eps]);
  const double BA_e_crf = BA_q * (BA_distisoft - BA_crf_2cut3i * dist2 - m_crf_cut[eps]);
  const double B_e_crf = B_q * (B_distisoft - B_crf_2cut3i * dist2 - m_crf_cut[eps]);

  e_crf = (lamA_crf * A_e_crf + lamAB_crf * AB_e_crf +
           lamBA_crf * BA_e_crf + lamB_crf * B_e_crf) * math::four_pi_eps_i;

  de_crf = ((A_e_crf * lamA_crf_deriv - (alpha_crf /2) * lamA_crf *
           (A_dist3isoft - A_crf_pert * dist2) * lamA_crfs_deriv ) * A_q +
           (AB_e_crf * lamAB_crf_deriv - (alpha_crf /2) * lamAB_crf *
           (AB_dist3isoft - AB_crf_pert * dist2) * lamAB_crfs_deriv ) * AB_q +
           (BA_e_crf * lamBA_crf_deriv - (alpha_crf /2) * lamBA_crf *
           (BA_dist3isoft - BA_crf_pert * dist2) * lamBA_crfs_deriv ) * BA_q  +
           (B_e_crf * lamB_crf_deriv - (alpha_crf /2) * lamB_crf *
           (B_dist3isoft - B_crf_pert * dist2) * lamB_crfs_deriv ) * B_q )
           * math::four_pi_eps_i;

  pert2d_lamI_crf = ((A_e_crf * -(1 - lamJ_crf)) * A_q + (AB_e_crf * -lamJ_crf ) * AB_q +
                  (BA_e_crf * (1 - lamJ_crf)) * BA_q + (B_e_crf * lamJ_crf) * B_q)
                   * math::four_pi_eps_i;

  pert2d_lamI_crfs = ( lamA_crf * A_q * (A_dist3isoft - A_crf_pert * dist2) * dfA_dlamI_crfs 
             + lamAB_crf * AB_q * (AB_dist3isoft - AB_crf_pert * dist2) * dfAB_dlamI_crfs
             + lamBA_crf * BA_q * (BA_dist3isoft - BA_crf_pert * dist2) * dfBA_dlamI_crfs
             + lamB_crf * B_q * (B_dist3isoft - B_crf_pert * dist2) * dfB_dlamI_crfs)
             * -(alpha_crf /2) * math::four_pi_eps_i;

  pert2d_lamJ_crf = ((A_e_crf * -(1 - lamI_crf)) * A_q + (AB_e_crf * (1 - lamI_crf)) * AB_q +
                  (BA_e_crf * -lamI_crf) * BA_q + (B_e_crf * lamI_crf) * B_q)
                   * math::four_pi_eps_i;

  pert2d_lamJ_crfs = ( lamA_crf * A_q * (A_dist3isoft - A_crf_pert * dist2) * dfA_dlamJ_crfs 
             + lamAB_crf * AB_q * (AB_dist3isoft - AB_crf_pert * dist2) * dfAB_dlamJ_crfs
             + lamBA_crf * BA_q * (BA_dist3isoft - BA_crf_pert * dist2) * dfBA_dlamJ_crfs
             + lamB_crf * B_q * (B_dist3isoft - B_crf_pert * dist2) * dfB_dlamJ_crfs)
             * -(alpha_crf /2) * math::four_pi_eps_i;

} //2D perturbation

inline void interaction::Perturbed_Nonbonded_Term
::lj_crf_soft_2d_other_interaction(math::Vec const &r,
        double const A_c6, double const A_c12,
        double const B_c6, double const B_c12,
        double const A_q, double const B_q,
        double const alpha_lj, double const alpha_crf,
        double &force1, double &force6, double &force12,
        double &e_lj, double &e_crf,
        double &de_lj, double & de_crf, 
        double &pert2d_lj, double &pert2d_ljs, double &pert2d_crf, double &pert2d_crfs,
        unsigned int eps) {

  double A_c126, B_c126;

  if (A_c6 != 0) A_c126 = A_c12 / A_c6;
  else A_c126 = 0.0;
  if (B_c6 != 0) B_c126 = B_c12 / B_c6;
  else B_c126 = 0.0;

  const double dist2 = abs2(r);
  assert(dist2 != 0);

  const double A_dist2soft = dist2 + alpha_crf*m_B_crfs_lambda2;
  const double B_dist2soft = dist2 + alpha_crf*m_A_crfs_lambda2;

  const double A_distisoft = 1.0 / sqrt(A_dist2soft);
  const double B_distisoft = 1.0 / sqrt(B_dist2soft);

  const double A_dist3isoft = A_distisoft / A_dist2soft;
  const double B_dist3isoft = B_distisoft / B_dist2soft;

  const double dist4 = dist2 * dist2;
  const double dist6 = dist4 * dist2;

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

  // substitute A_dist3isoft thing. just like here -- daniel
  force1 = (m_A_crf_lambda_n * A_q * (A_dist3isoft + A_crf_cut3i) +
          m_B_crf_lambda_n * B_q * (B_dist3isoft + B_crf_cut3i)) *
          math::four_pi_eps_i;

  force6 = -6.0 * (m_A_lj_lambda_n * A_c6 * A_dist6isoft * A_dist6isoft +
          m_B_lj_lambda_n * B_c6 * B_dist6isoft * B_dist6isoft) * dist4;

  force12 = 12 * (m_A_lj_lambda_n * A_c12 * A_dist6isoft * A_dist6isoft * A_dist6isoft +
          m_B_lj_lambda_n * B_c12 * B_dist6isoft * B_dist6isoft * B_dist6isoft) * dist4;

  const double A_e_lj = (A_c12 * A_dist6isoft - A_c6) * A_dist6isoft;
  const double B_e_lj = (B_c12 * B_dist6isoft - B_c6) * B_dist6isoft;

  const double A_e_crf = A_q * (A_distisoft - A_crf_2cut3i * dist2 - m_crf_cut[eps]);
  const double B_e_crf = B_q * (B_distisoft - B_crf_2cut3i * dist2 - m_crf_cut[eps]);

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

//  std::cout << "\tA_e_lj " <<  A_e_lj << " m_A_lj_lambda_n " <<  m_A_lj_lambda_n <<
//               " B_e_lj " << B_e_lj << " m_B_lj_lambda_n " << m_B_lj_lambda_n << std::endl;

  e_lj = m_A_lj_lambda_n * A_e_lj + m_B_lj_lambda_n * B_e_lj;

  e_crf = (m_A_crf_lambda_n * A_e_crf + m_B_crf_lambda_n * B_e_crf) * math::four_pi_eps_i;

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

  pert2d_lj = m_lambda_exp * (pow(lj_lambda, m_lambda_exp-1) * B_e_lj - pow(1-lj_lambda, m_lambda_exp-1) * A_e_lj);
  pert2d_ljs = -2.0 * alpha_lj * (m_A_lj_lambda_n * ljs_lambda * A_c126 * A_dist6isoft * A_dist6isoft *
          (2 * A_c12 * A_dist6isoft - A_c6) -
          m_B_lj_lambda_n * (1-ljs_lambda) * B_c126 * B_dist6isoft * B_dist6isoft *
          (2 * B_c12 * B_dist6isoft - B_c6));

  pert2d_crf = m_lambda_exp * (pow(crf_lambda,m_lambda_exp-1) * B_e_crf - pow(1-crf_lambda,m_lambda_exp-1)* A_e_crf) *
          math::four_pi_eps_i;
  pert2d_crfs = -(m_A_crf_lambda_n * A_q * crfs_lambda * (A_dist3isoft - A_crf_pert * dist2) -
          m_B_crf_lambda_n * B_q * (1-crfs_lambda) * (B_dist3isoft - B_crf_pert * dist2)) *
          math::four_pi_eps_i * alpha_crf;
}


inline void interaction::Perturbed_Nonbonded_Term
::pol_lj_crf_soft_interaction(math::Vec const &r, math::Vec const &rp1,
        math::Vec const &rp2, math::Vec const &rpp,
        double const A_c6, double const A_c12,
        double const B_c6, double const B_c12,
        double const A_qi, double const B_qi,
        double const A_qj, double const B_qj,
        double const cqi, double const cqj,
        double const alpha_lj, double const alpha_crf,
        double force1[],
        double &force6, double &force12,
        double &e_lj, double &e_crf,
        double &de_lj, double & de_crf, unsigned int eps) {
  double A_c126, B_c126;

  if (A_c6 != 0) A_c126 = A_c12 / A_c6;
  else A_c126 = 0.0;
  if (B_c6 != 0) B_c126 = B_c12 / B_c6;
  else B_c126 = 0.0;

  const double A_qeps = (A_qi - cqi) * (A_qj - cqj) * math::four_pi_eps_i;
  const double A_qepsp1 = (A_qi - cqi) * cqj * math::four_pi_eps_i;
  const double A_qepsp2 = (A_qj - cqj) * cqi * math::four_pi_eps_i;
  const double A_qepspp = cqi * cqj * math::four_pi_eps_i;

  const double B_qeps = (B_qi - cqi) * (B_qj - cqj) * math::four_pi_eps_i;
  const double B_qepsp1 = (B_qi - cqi) * cqj * math::four_pi_eps_i;
  const double B_qepsp2 = (B_qj - cqj) * cqi * math::four_pi_eps_i;
  const double B_qepspp = cqi * cqj * math::four_pi_eps_i;

  const double dist2 = abs2(r);
  const double dist2p1 = abs2(rp1);
  const double dist2p2 = abs2(rp2);
  const double dist2pp = abs2(rpp);
  assert(dist2 != 0);

  const double A_al2 = alpha_crf*m_B_crfs_lambda2;
  const double B_al2 = alpha_crf*m_A_crfs_lambda2;
  const double A_dist2soft = dist2 + A_al2;
  const double A_dist2p1soft = dist2p1 + A_al2;
  const double A_dist2p2soft = dist2p2 + A_al2;
  const double A_dist2ppsoft = dist2pp + A_al2;
  const double B_dist2soft = dist2 + B_al2;
  const double B_dist2p1soft = dist2p1 + B_al2;
  const double B_dist2p2soft = dist2p2 + B_al2;
  const double B_dist2ppsoft = dist2pp + B_al2;

  const double A_distisoft = 1.0 / sqrt(A_dist2soft);
  const double A_distip1soft = 1.0 / sqrt(A_dist2p1soft);
  const double A_distip2soft = 1.0 / sqrt(A_dist2p2soft);
  const double A_distippsoft = 1.0 / sqrt(A_dist2ppsoft);
  const double B_distisoft = 1.0 / sqrt(B_dist2soft);
  const double B_distip1soft = 1.0 / sqrt(B_dist2p1soft);
  const double B_distip2soft = 1.0 / sqrt(B_dist2p2soft);
  const double B_distippsoft = 1.0 / sqrt(B_dist2ppsoft);

  const double A_dist3isoft = A_distisoft / A_dist2soft;
  const double A_dist3ip1soft = A_distip1soft / A_dist2p1soft;
  const double A_dist3ip2soft = A_distip2soft / A_dist2p2soft;
  const double A_dist3ippsoft = A_distippsoft / A_dist2ppsoft;
  const double B_dist3isoft = B_distisoft / B_dist2soft;
  const double B_dist3ip1soft = B_distip1soft / B_dist2p1soft;
  const double B_dist3ip2soft = B_distip2soft / B_dist2p2soft;
  const double B_dist3ippsoft = B_distippsoft / B_dist2ppsoft;

  const double dist4 = dist2 * dist2;
  const double dist6 = dist4 * dist2;

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

  // substitute A_dist3isoft thing. just like here -- daniel
  force1[0] = m_A_crf_lambda_n * A_qeps * (A_dist3isoft + A_crf_cut3i) +
          m_B_crf_lambda_n * B_qeps * (B_dist3isoft + B_crf_cut3i);
  force1[1] = m_A_crf_lambda_n * A_qepsp1 * (A_dist3ip1soft + A_crf_cut3i) +
          m_B_crf_lambda_n * B_qepsp1 * (B_dist3ip1soft + B_crf_cut3i);
  force1[2] = m_A_crf_lambda_n * A_qepsp2 * (A_dist3ip2soft + A_crf_cut3i) +
          m_B_crf_lambda_n * B_qepsp2 * (B_dist3ip2soft + B_crf_cut3i);
  force1[3] = m_A_crf_lambda_n * A_qepspp * (A_dist3ippsoft + A_crf_cut3i) +
          m_B_crf_lambda_n * B_qepspp * (B_dist3ippsoft + B_crf_cut3i);

  force6 = -6.0 * (m_A_lj_lambda_n * A_c6 * A_dist6isoft * A_dist6isoft +
          m_B_lj_lambda_n * B_c6 * B_dist6isoft * B_dist6isoft) * dist4;

  force12 = 12 * (m_A_lj_lambda_n * A_c12 * A_dist6isoft * A_dist6isoft * A_dist6isoft +
          m_B_lj_lambda_n * B_c12 * B_dist6isoft * B_dist6isoft * B_dist6isoft) * dist4;

  const double A_e_lj = (A_c12 * A_dist6isoft - A_c6) * A_dist6isoft;
  const double B_e_lj = (B_c12 * B_dist6isoft - B_c6) * B_dist6isoft;

  double A_ecrf0 = A_qeps * (A_distisoft - A_crf_2cut3i * dist2 - m_crf_cut[eps]);
  double A_ecrf1 = A_qepsp1 * (A_distip1soft - A_crf_2cut3i * dist2p1 - m_crf_cut[eps]);
  double A_ecrf2 = A_qepsp2 * (A_distip2soft - A_crf_2cut3i * dist2p2 - m_crf_cut[eps]);
  double A_ecrfp = A_qepspp * (A_distippsoft - A_crf_2cut3i * dist2pp - m_crf_cut[eps]);
  double B_ecrf0 = B_qeps * (B_distisoft - B_crf_2cut3i * dist2 - m_crf_cut[eps]);
  double B_ecrf1 = B_qepsp1 * (B_distip1soft - B_crf_2cut3i * dist2p1 - m_crf_cut[eps]);
  double B_ecrf2 = B_qepsp2 * (B_distip2soft - B_crf_2cut3i * dist2p2 - m_crf_cut[eps]);
  double B_ecrfp = B_qepspp * (B_distippsoft - B_crf_2cut3i * dist2pp - m_crf_cut[eps]);

  const double A_e_crf = A_ecrf0 + A_ecrf1 + A_ecrf2 + A_ecrfp;
  const double B_e_crf = B_ecrf0 + B_ecrf1 + B_ecrf2 + B_ecrfp;

  e_lj = m_A_lj_lambda_n * A_e_lj + m_B_lj_lambda_n * B_e_lj;

  e_crf = (m_A_crf_lambda_n * A_e_crf + m_B_crf_lambda_n * B_e_crf);

  de_lj = -2.0 * alpha_lj * (m_A_lj_lambda_n * m_B_ljs_lambda * A_c126 * A_dist6isoft * A_dist6isoft * (2 * A_c12 * A_dist6isoft - A_c6) -
          m_B_lj_lambda_n * m_A_ljs_lambda * B_c126 * B_dist6isoft * B_dist6isoft * (2 * B_c12 * B_dist6isoft - B_c6))
          + m_lambda_exp * (m_B_lj_lambda_n_1 * B_e_lj - m_A_lj_lambda_n_1 * A_e_lj);

  de_crf = (-(m_A_crf_lambda_n * A_qeps * m_B_crfs_lambda * (A_dist3isoft - A_crf_pert * dist2) -
          m_B_crf_lambda_n * B_qeps * m_A_crfs_lambda * (B_dist3isoft - B_crf_pert * dist2)) * alpha_crf
          + m_lambda_exp * (m_B_crf_lambda_n_1 * B_ecrf0 - m_A_crf_lambda_n_1 * A_ecrf0)) +
          (-(m_A_crf_lambda_n * A_qepsp1 * m_B_crfs_lambda * (A_dist3ip1soft
          - A_crf_pert * dist2p1) - m_B_crf_lambda_n * B_qepsp1 * m_A_crfs_lambda *
          (B_dist3ip1soft - B_crf_pert * dist2p1)) * alpha_crf
          + m_lambda_exp * (m_B_crf_lambda_n_1 * B_ecrf1 - m_A_crf_lambda_n_1 * A_ecrf1)) +
          (-(m_A_crf_lambda_n * A_qepsp2 * m_B_crfs_lambda * (A_dist3ip2soft
          - A_crf_pert * dist2p2) - m_B_crf_lambda_n * B_qepsp2 * m_A_crfs_lambda *
          (B_dist3ip2soft - B_crf_pert * dist2p2)) * alpha_crf
          + m_lambda_exp * (m_B_crf_lambda_n_1 * B_ecrf2 - m_A_crf_lambda_n_1 * A_ecrf2)) +
          (-(m_A_crf_lambda_n * A_qepspp * m_B_crfs_lambda * (A_dist3ippsoft
          - A_crf_pert * dist2pp) - m_B_crf_lambda_n * B_qepspp * m_A_crfs_lambda *
          (B_dist3ippsoft - B_crf_pert * dist2pp)) * alpha_crf
          + m_lambda_exp * (m_B_crf_lambda_n_1 * B_ecrfp - m_A_crf_lambda_n_1 * A_ecrfp));
}

inline void interaction::Perturbed_Nonbonded_Term
::pol_off_lj_crf_soft_interaction(math::Vec const &r, math::Vec const &rm, math::Vec const &rp1,
        math::Vec const &rp2, math::Vec const &rpp,
        double const A_c6, double const A_c12,
        double const B_c6, double const B_c12,
        double const A_qi, double const B_qi,
        double const A_qj, double const B_qj,
        double const cqi, double const cqj,
        double const alpha_lj, double const alpha_crf,
        double force1[],
        double &force6, double &force12,
        double &e_lj, double &e_crf,
        double &de_lj, double & de_crf, unsigned int eps) {
  double A_c126, B_c126;

  if (A_c6 != 0) A_c126 = A_c12 / A_c6;
  else A_c126 = 0.0;
  if (B_c6 != 0) B_c126 = B_c12 / B_c6;
  else B_c126 = 0.0;

  const double A_qeps = (A_qi - cqi) * (A_qj - cqj) * math::four_pi_eps_i;
  const double A_qepsp1 = (A_qi - cqi) * cqj * math::four_pi_eps_i;
  const double A_qepsp2 = (A_qj - cqj) * cqi * math::four_pi_eps_i;
  const double A_qepspp = cqi * cqj * math::four_pi_eps_i;

  const double B_qeps = (B_qi - cqi) * (B_qj - cqj) * math::four_pi_eps_i;
  const double B_qepsp1 = (B_qi - cqi) * cqj * math::four_pi_eps_i;
  const double B_qepsp2 = (B_qj - cqj) * cqi * math::four_pi_eps_i;
  const double B_qepspp = cqi * cqj * math::four_pi_eps_i;

  const double dist2 = abs2(r);
  const double dist2m = abs2(rm);
  const double dist2p1 = abs2(rp1);
  const double dist2p2 = abs2(rp2);
  const double dist2pp = abs2(rpp);
  assert(dist2 != 0);

  const double A_al2 = alpha_crf*m_B_crfs_lambda2;
  const double B_al2 = alpha_crf*m_A_crfs_lambda2;
  const double A_dist2soft = dist2m + A_al2;
  const double A_dist2p1soft = dist2p1 + A_al2;
  const double A_dist2p2soft = dist2p2 + A_al2;
  const double A_dist2ppsoft = dist2pp + A_al2;
  const double B_dist2soft = dist2m + B_al2;
  const double B_dist2p1soft = dist2p1 + B_al2;
  const double B_dist2p2soft = dist2p2 + B_al2;
  const double B_dist2ppsoft = dist2pp + B_al2;

  const double A_distisoft = 1.0 / sqrt(A_dist2soft);
  const double A_distip1soft = 1.0 / sqrt(A_dist2p1soft);
  const double A_distip2soft = 1.0 / sqrt(A_dist2p2soft);
  const double A_distippsoft = 1.0 / sqrt(A_dist2ppsoft);
  const double B_distisoft = 1.0 / sqrt(B_dist2soft);
  const double B_distip1soft = 1.0 / sqrt(B_dist2p1soft);
  const double B_distip2soft = 1.0 / sqrt(B_dist2p2soft);
  const double B_distippsoft = 1.0 / sqrt(B_dist2ppsoft);

  const double A_dist3isoft = A_distisoft / A_dist2soft;
  const double A_dist3ip1soft = A_distip1soft / A_dist2p1soft;
  const double A_dist3ip2soft = A_distip2soft / A_dist2p2soft;
  const double A_dist3ippsoft = A_distippsoft / A_dist2ppsoft;
  const double B_dist3isoft = B_distisoft / B_dist2soft;
  const double B_dist3ip1soft = B_distip1soft / B_dist2p1soft;
  const double B_dist3ip2soft = B_distip2soft / B_dist2p2soft;
  const double B_dist3ippsoft = B_distippsoft / B_dist2ppsoft;

  const double dist4 = dist2 * dist2;
  const double dist6 = dist4 * dist2;

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

  // substitute A_dist3isoft thing. just like here -- daniel
  force1[0] = m_A_crf_lambda_n * A_qeps * (A_dist3isoft + A_crf_cut3i) +
          m_B_crf_lambda_n * B_qeps * (B_dist3isoft + B_crf_cut3i);
  force1[1] = m_A_crf_lambda_n * A_qepsp1 * (A_dist3ip1soft + A_crf_cut3i) +
          m_B_crf_lambda_n * B_qepsp1 * (B_dist3ip1soft + B_crf_cut3i);
  force1[2] = m_A_crf_lambda_n * A_qepsp2 * (A_dist3ip2soft + A_crf_cut3i) +
          m_B_crf_lambda_n * B_qepsp2 * (B_dist3ip2soft + B_crf_cut3i);
  force1[3] = m_A_crf_lambda_n * A_qepspp * (A_dist3ippsoft + A_crf_cut3i) +
          m_B_crf_lambda_n * B_qepspp * (B_dist3ippsoft + B_crf_cut3i);

  force6 = -6.0 * (m_A_lj_lambda_n * A_c6 * A_dist6isoft * A_dist6isoft +
          m_B_lj_lambda_n * B_c6 * B_dist6isoft * B_dist6isoft) * dist4;

  force12 = 12 * (m_A_lj_lambda_n * A_c12 * A_dist6isoft * A_dist6isoft * A_dist6isoft +
          m_B_lj_lambda_n * B_c12 * B_dist6isoft * B_dist6isoft * B_dist6isoft) * dist4;

  const double A_e_lj = (A_c12 * A_dist6isoft - A_c6) * A_dist6isoft;
  const double B_e_lj = (B_c12 * B_dist6isoft - B_c6) * B_dist6isoft;

  double A_ecrf0 = A_qeps * (A_distisoft - A_crf_2cut3i * dist2m - m_crf_cut[eps]);
  double A_ecrf1 = A_qepsp1 * (A_distip1soft - A_crf_2cut3i * dist2p1 - m_crf_cut[eps]);
  double A_ecrf2 = A_qepsp2 * (A_distip2soft - A_crf_2cut3i * dist2p2 - m_crf_cut[eps]);
  double A_ecrfp = A_qepspp * (A_distippsoft - A_crf_2cut3i * dist2pp - m_crf_cut[eps]);
  double B_ecrf0 = B_qeps * (B_distisoft - B_crf_2cut3i * dist2m - m_crf_cut[eps]);
  double B_ecrf1 = B_qepsp1 * (B_distip1soft - B_crf_2cut3i * dist2p1 - m_crf_cut[eps]);
  double B_ecrf2 = B_qepsp2 * (B_distip2soft - B_crf_2cut3i * dist2p2 - m_crf_cut[eps]);
  double B_ecrfp = B_qepspp * (B_distippsoft - B_crf_2cut3i * dist2pp - m_crf_cut[eps]);

  const double A_e_crf = A_ecrf0 + A_ecrf1 + A_ecrf2 + A_ecrfp;
  const double B_e_crf = B_ecrf0 + B_ecrf1 + B_ecrf2 + B_ecrfp;

  e_lj = m_A_lj_lambda_n * A_e_lj + m_B_lj_lambda_n * B_e_lj;

  e_crf = (m_A_crf_lambda_n * A_e_crf + m_B_crf_lambda_n * B_e_crf);

  de_lj = -2.0 * alpha_lj * (m_A_lj_lambda_n * m_B_ljs_lambda * A_c126 * A_dist6isoft * A_dist6isoft *
          (2 * A_c12 * A_dist6isoft - A_c6) -
          m_B_lj_lambda_n * m_A_ljs_lambda * B_c126 * B_dist6isoft * B_dist6isoft *
          (2 * B_c12 * B_dist6isoft - B_c6))
          + m_lambda_exp * (m_B_lj_lambda_n_1 * B_e_lj - m_A_lj_lambda_n_1 * A_e_lj);

  de_crf = (-(m_A_crf_lambda_n * A_qeps * m_B_crfs_lambda * (A_dist3isoft - A_crf_pert * dist2m) -
          m_B_crf_lambda_n * B_qeps * m_A_crfs_lambda * (B_dist3isoft - B_crf_pert * dist2m)) * alpha_crf
          + m_lambda_exp * (m_B_crf_lambda_n_1 * B_ecrf0 - m_A_crf_lambda_n_1 * A_ecrf0)) +
          (-(m_A_crf_lambda_n * A_qepsp1 * m_B_crfs_lambda * (A_dist3ip1soft
          - A_crf_pert * dist2p1) - m_B_crf_lambda_n * B_qepsp1 * m_A_crfs_lambda *
          (B_dist3ip1soft - B_crf_pert * dist2p1)) * alpha_crf
          + m_lambda_exp * (m_B_crf_lambda_n_1 * B_ecrf1 - m_A_crf_lambda_n_1 * A_ecrf1)) +
          (-(m_A_crf_lambda_n * A_qepsp2 * m_B_crfs_lambda * (A_dist3ip2soft
          - A_crf_pert * dist2p2) - m_B_crf_lambda_n * B_qepsp2 * m_A_crfs_lambda *
          (B_dist3ip2soft - B_crf_pert * dist2p2)) * alpha_crf
          + m_lambda_exp * (m_B_crf_lambda_n_1 * B_ecrf2 - m_A_crf_lambda_n_1 * A_ecrf2)) +
          (-(m_A_crf_lambda_n * A_qepspp * m_B_crfs_lambda * (A_dist3ippsoft
          - A_crf_pert * dist2pp) - m_B_crf_lambda_n * B_qepspp * m_A_crfs_lambda *
          (B_dist3ippsoft - B_crf_pert * dist2pp)) * alpha_crf
          + m_lambda_exp * (m_B_crf_lambda_n_1 * B_ecrfp - m_A_crf_lambda_n_1 * A_ecrfp));
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
        double &de_lj, double & de_crf, unsigned int eps) {
  double A_c126, B_c126;

  if (A_c6 != 0) A_c126 = A_c12 / A_c6;
  else A_c126 = 0.0;
  if (B_c6 != 0) B_c126 = B_c12 / B_c6;
  else B_c126 = 0.0;

  const double dist2 = abs2(r);
  assert(dist2 != 0);

  const double A_dist2soft = dist2 + alpha_crf*m_B_crfs_lambda2;
  const double B_dist2soft = dist2 + alpha_crf*m_A_crfs_lambda2;

  const double A_distisoft = 1.0 / sqrt(A_dist2soft);
  const double B_distisoft = 1.0 / sqrt(B_dist2soft);

  const double A_dist3isoft = A_distisoft / A_dist2soft;
  const double B_dist3isoft = B_distisoft / B_dist2soft;

  const double dist4 = dist2 * dist2;
  const double dist6 = dist4 * dist2;

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

  const double A_lj_scaled_lambda_n = A_scale * m_A_lj_lambda_n;
  const double B_lj_scaled_lambda_n = B_scale * m_B_lj_lambda_n;
  const double A_crf_scaled_lambda_n = A_scale * m_A_crf_lambda_n;
  const double B_crf_scaled_lambda_n = B_scale * m_B_crf_lambda_n;


  force1 = (A_crf_scaled_lambda_n * A_q * (A_dist3isoft + A_crf_cut3i) +
          B_crf_scaled_lambda_n * B_q * (B_dist3isoft + B_crf_cut3i)) *
          math::four_pi_eps_i;

  force6 = -6.0 * (A_lj_scaled_lambda_n * A_c6 * A_dist6isoft * A_dist6isoft +
          B_lj_scaled_lambda_n * B_c6 * B_dist6isoft * B_dist6isoft) * dist4;

  force12 = 12 * (A_lj_scaled_lambda_n * A_c12 * A_dist6isoft * A_dist6isoft * A_dist6isoft +
          B_lj_scaled_lambda_n * B_c12 * B_dist6isoft * B_dist6isoft * B_dist6isoft) * dist4;

  const double A_e_lj = (A_c12 * A_dist6isoft - A_c6) * A_dist6isoft;
  const double B_e_lj = (B_c12 * B_dist6isoft - B_c6) * B_dist6isoft;

  const double A_e_crf = A_q * (A_distisoft - A_crf_2cut3i * dist2 - m_crf_cut[eps]);
  const double B_e_crf = B_q * (B_distisoft - B_crf_2cut3i * dist2 - m_crf_cut[eps]);


  e_lj = A_lj_scaled_lambda_n * A_e_lj + B_lj_scaled_lambda_n * B_e_lj;

  e_crf = (A_crf_scaled_lambda_n * A_e_crf + B_crf_scaled_lambda_n * B_e_crf) *
          math::four_pi_eps_i;

  de_lj = -2.0 * alpha_lj * (A_lj_scaled_lambda_n * m_B_ljs_lambda * A_c126 *
          A_dist6isoft * A_dist6isoft * (2 * A_c12 * A_dist6isoft - A_c6) -
          B_lj_scaled_lambda_n * m_A_ljs_lambda * B_c126 *
          B_dist6isoft * B_dist6isoft * (2 * B_c12 * B_dist6isoft - B_c6))
          + m_lambda_exp * (B_scale * m_B_lj_lambda_n_1 * B_e_lj - A_scale * m_A_lj_lambda_n_1 * A_e_lj);

  de_crf = (-alpha_crf * (A_crf_scaled_lambda_n * A_q * m_B_crfs_lambda *
          (A_dist3isoft - A_crf_pert * dist2) -
          B_crf_scaled_lambda_n * B_q * m_A_crfs_lambda *
          (B_dist3isoft - B_crf_pert * dist2))
          + m_lambda_exp * (B_scale * m_B_crf_lambda_n_1 * B_e_crf -
          A_scale * m_A_crf_lambda_n_1 * A_e_crf))
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
        double &de_lj, double & de_crf) {
  double A_c126, B_c126;

  if (A_c6 != 0) A_c126 = A_c12 / A_c6;
  else A_c126 = 0.0;
  if (B_c6 != 0) B_c126 = B_c12 / B_c6;
  else B_c126 = 0.0;

  const double dist2 = abs2(r);
  assert(dist2 != 0);

  const double dist = sqrt(dist2);

  const double A_dist2soft = dist2 + alpha_crf*m_B_crfs_lambda2;
  const double B_dist2soft = dist2 + alpha_crf*m_A_crfs_lambda2;

  const double A_distisoft = 1.0 / sqrt(A_dist2soft);
  const double B_distisoft = 1.0 / sqrt(B_dist2soft);

  const double A_dist3isoft = A_distisoft / A_dist2soft;
  const double B_dist3isoft = B_distisoft / B_dist2soft;

  const double A_dist2soft_cut = nb_cutoff * nb_cutoff + alpha_crf * m_B_crfs_lambda2;
  const double B_dist2soft_cut = nb_cutoff * nb_cutoff + alpha_crf * m_A_crfs_lambda2;

  const double A_distisoft_cut = 1.0 / sqrt(A_dist2soft_cut);
  const double B_distisoft_cut = 1.0 / sqrt(B_dist2soft_cut);

  const double A_dist3isoft_cut = A_distisoft_cut / A_dist2soft_cut;
  const double B_dist3isoft_cut = B_distisoft_cut / B_dist2soft_cut;

  const double dist4 = dist2 * dist2;
  const double dist6 = dist4 * dist2;

  const double A_dist6soft = dist6 + alpha_lj * m_B_ljs_lambda2*A_c126;
  const double B_dist6soft = dist6 + alpha_lj * m_A_ljs_lambda2*B_c126;

  const double A_dist6isoft = 1.0 / A_dist6soft;
  const double B_dist6isoft = 1.0 / B_dist6soft;

  // setting the cgrain C constant here
  // for the cgrain, the const C depends on lambda

  const double A_distsoft_cut = (pow(nb_cutoff, 6)
          + A_c126 * alpha_lj * m_B_ljs_lambda2);
  const double B_distsoft_cut = (pow(nb_cutoff, 6)
          + B_c126 * alpha_lj * m_A_ljs_lambda2);

  const double A_dist6isoft_cut = 1.0 / A_distsoft_cut;
  const double B_dist6isoft_cut = 1.0 / B_distsoft_cut;

  // state A
  A_C_cg12 = A_dist6isoft_cut * A_dist6isoft_cut
          - A_cg12 / 3 * nb_cutoff * nb_cutoff * nb_cutoff
          - B_cg12 / 4 * nb_cutoff * nb_cutoff * nb_cutoff * nb_cutoff;
  A_C_cg6 = A_dist6isoft_cut
          - A_cg6 / 3 * nb_cutoff * nb_cutoff * nb_cutoff
          - B_cg6 / 4 * nb_cutoff * nb_cutoff * nb_cutoff * nb_cutoff;
  A_C_cg1 = 1.0 / sqrt((nb_cutoff * nb_cutoff + alpha_crf * m_B_crfs_lambda2))
          - A_cg1 / 3 * nb_cutoff * nb_cutoff * nb_cutoff
          - B_cg1 / 4 * nb_cutoff * nb_cutoff * nb_cutoff * nb_cutoff;

  // state B
  B_C_cg12 = B_dist6isoft_cut * B_dist6isoft_cut
          - A_cg12 / 3 * nb_cutoff * nb_cutoff * nb_cutoff
          - B_cg12 / 4 * nb_cutoff * nb_cutoff * nb_cutoff * nb_cutoff;
  B_C_cg6 = B_dist6isoft_cut
          - A_cg6 / 3 * nb_cutoff * nb_cutoff * nb_cutoff
          - B_cg6 / 4 * nb_cutoff * nb_cutoff * nb_cutoff * nb_cutoff;
  B_C_cg1 = 1.0 / sqrt((nb_cutoff * nb_cutoff + alpha_crf * m_A_crfs_lambda2))
          - A_cg1 / 3 * nb_cutoff * nb_cutoff * nb_cutoff
          - B_cg1 / 4 * nb_cutoff * nb_cutoff * nb_cutoff * nb_cutoff;


  DEBUG(11, "just checking\nm_A_lambda " << m_A_ljs_lambda
          << "\nm_B_lambda " << m_B_ljs_lambda
          << "\nm_A_lambda2 " << m_A_ljs_lambda2
          << "\nm_B_lambda2 " << m_B_ljs_lambda2
          << "\nm_A_lambda_n " << m_A_lj_lambda_n
          << "\nm_B_lambda_n" << m_B_lj_lambda_n
          << "\nm_A_lambda_n_1 " << m_A_lj_lambda_n_1
          << "\nm_B_lambda_n_1 " << m_B_lj_lambda_n_1);

  // coulomb 
  force1 = (m_A_crf_lambda_n * A_q * (A_dist3isoft +
          A_cg1 * dist +
          B_cg1 * dist2) +
          m_B_crf_lambda_n * B_q * (B_dist3isoft +
          A_cg1 * dist +
          B_cg1 * dist2)) *
          math::four_pi_eps_i / cgrain_eps[0];

  const double A_e_crf = A_q * (A_distisoft
          - A_cg1 / 3 * dist2 * dist
          - B_cg1 / 4 * dist2 * dist2
          - A_C_cg1);
  const double B_e_crf = B_q * (B_distisoft
          - A_cg1 / 3 * dist2 * dist
          - B_cg1 / 4 * dist2 * dist2
          - B_C_cg1);

  e_crf = (m_A_crf_lambda_n * A_e_crf + m_B_crf_lambda_n * B_e_crf) * math::four_pi_eps_i / cgrain_eps[0];

  de_crf = -(m_A_crf_lambda_n * A_q * m_B_crfs_lambda * (A_dist3isoft - A_dist3isoft_cut) -
          m_B_crf_lambda_n * B_q * m_A_crfs_lambda * (B_dist3isoft - B_dist3isoft_cut))
          * math::four_pi_eps_i * alpha_crf / cgrain_eps[0]
          + m_lambda_exp * (m_B_crf_lambda_n_1 * B_e_crf - m_A_crf_lambda_n_1 * A_e_crf)
          * math::four_pi_eps_i / cgrain_eps[0];

  //LJ
  // after some consideration and an appropriate measure of doubt we
  // change the second 6.0 * into a 1.0 *
  force6 = -6.0 * (m_A_lj_lambda_n * A_c6 * A_dist6isoft * A_dist6isoft +
          m_B_lj_lambda_n * B_c6 * B_dist6isoft * B_dist6isoft) * dist4
          - 1.0 * (m_A_lj_lambda_n * A_c6 * (A_cg6 * dist2 + B_cg6 * dist2 * dist) +
          m_B_lj_lambda_n * B_c6 * (A_cg6 * dist2 + B_cg6 * dist2 * dist)) / dist;

  force12 = 12.0 * (m_A_lj_lambda_n * A_c12 * A_dist6isoft * A_dist6isoft * A_dist6isoft +
          m_B_lj_lambda_n * B_c12 * B_dist6isoft * B_dist6isoft * B_dist6isoft) * dist4 +
          +1.0 * (m_A_lj_lambda_n * A_c12 * (A_cg12 * dist2 + B_cg12 * dist2 * dist) +
          m_B_lj_lambda_n * B_c12 * (A_cg12 * dist2 + B_cg12 * dist2 * dist)) / dist;


  const double A_e_lj = A_c12 * (A_dist6isoft * A_dist6isoft
          - A_cg12 / 3 * dist2 * dist
          - B_cg12 / 4 * dist2 * dist2
          - A_C_cg12)
          - A_c6 * (A_dist6isoft
          - A_cg6 / 3 * dist2 * dist
          - B_cg6 / 4 * dist2 * dist2
          - A_C_cg6);

  const double B_e_lj = B_c12 * (B_dist6isoft * B_dist6isoft
          - A_cg12 / 3 * dist2 * dist
          - B_cg12 / 4 * dist2 * dist2
          - B_C_cg12)
          - B_c6 * (B_dist6isoft
          - A_cg6 / 3 * dist2 * dist
          - B_cg6 / 4 * dist2 * dist2
          - B_C_cg6);

  e_lj = m_A_lj_lambda_n * A_e_lj + m_B_lj_lambda_n * B_e_lj;

  // there is a bug here from the previous version
  // the constant C also depends on lambda (correct for coulomb)
  // we have then to subtract the term at the cut off distance
  de_lj = -2.0 * alpha_lj * (m_A_lj_lambda_n * m_B_ljs_lambda * A_c126 * A_dist6isoft * A_dist6isoft *
          (2 * A_c12 * A_dist6isoft - A_c6) -
          m_B_lj_lambda_n * m_A_ljs_lambda * B_c126 * B_dist6isoft * B_dist6isoft *
          (2 * B_c12 * B_dist6isoft - B_c6))
          + 2.0 * alpha_lj * (m_A_lj_lambda_n * m_B_ljs_lambda * A_c126 *
          A_dist6isoft_cut * A_dist6isoft_cut *
          (2 * A_c12 * A_dist6isoft_cut - A_c6) -
          m_B_lj_lambda_n * m_A_ljs_lambda * B_c126 *
          B_dist6isoft_cut * B_dist6isoft_cut *
          (2 * B_c12 * B_dist6isoft_cut - B_c6))
          + m_lambda_exp * (m_B_lj_lambda_n_1 * B_e_lj - m_A_lj_lambda_n_1 * A_e_lj);

  DEBUG(11, "cgrain_pert\nr_ij " << dist
          << "\nf1 " << force1
          << "\ne_c " << e_crf
          << "\nde_c " << de_crf
          << "\nf6 " << force6
          << "\nf12 " << force12
          << "\nf12_6" << force12 + force6
          << "\ne_lj " << e_lj
          << "\nde_lj " << de_lj);
}

/**
 * helper function to calculate a term of the electric field 
 * at a given position for the polarisation
 */
inline void interaction::Perturbed_Nonbonded_Term
::electric_field_soft_interaction(math::Vec const &r,
        math::Vec const &rprime,
        double const alpha_crf,
        double A_qj, double B_qj, double cgj,
        math::Vec &e_el, unsigned int eps) {

  DEBUG(14, "\t\tenergy field term for polarisation");

  math::Vec A_el, B_el;

  assert(abs2(r) != 0);
  const double dist2j = abs2(r);
  const double dist2p = abs2(rprime);

  const double A_dist2soft = dist2j + alpha_crf*m_B_crfs_lambda2;
  const double A_dist2psoft = dist2p + alpha_crf*m_B_crfs_lambda2;
  const double B_dist2soft = dist2j + alpha_crf*m_A_crfs_lambda2;
  const double B_dist2psoft = dist2p + alpha_crf*m_A_crfs_lambda2;

  const double A_dist3isoft = 1.0 / sqrt(A_dist2soft * A_dist2soft * A_dist2soft);
  const double A_dist3ipsoft = 1.0 / sqrt(A_dist2psoft * A_dist2psoft * A_dist2psoft);
  const double B_dist3isoft = 1.0 / sqrt(B_dist2soft * B_dist2soft * B_dist2soft);
  const double B_dist3ipsoft = 1.0 / sqrt(B_dist2psoft * B_dist2psoft * B_dist2psoft);

  const double A_qeps = (A_qj - cgj) * math::four_pi_eps_i;
  const double B_qeps = (B_qj - cgj) * math::four_pi_eps_i;
  const double qepsp = cgj * math::four_pi_eps_i;

  const double A_cut2soft = m_cut2 + alpha_crf * m_B_crfs_lambda2;
  const double B_cut2soft = m_cut2 + alpha_crf * m_A_crfs_lambda2;

  const double A_cut2soft3 = A_cut2soft * A_cut2soft * A_cut2soft;
  const double B_cut2soft3 = B_cut2soft * B_cut2soft * B_cut2soft;

  const double A_crf_cut3i = m_crf[eps] / sqrt(A_cut2soft3);
  const double B_crf_cut3i = m_crf[eps] / sqrt(B_cut2soft3);

  A_el = A_qeps * (A_dist3isoft + A_crf_cut3i) * r + qepsp * (A_dist3ipsoft + A_crf_cut3i) * rprime;
  B_el = B_qeps * (B_dist3isoft + B_crf_cut3i) * r + qepsp * (B_dist3ipsoft + B_crf_cut3i) * rprime;

  e_el = m_A_crf_lambda_n * A_el + m_B_crf_lambda_n*B_el;
}

/**
 * helper function to calculate the self energy 
 * at a given atom.
 */
inline void interaction::Perturbed_Nonbonded_Term
::self_energy_soft_interaction(double A_alpha, double B_alpha,
        double e_i2, double &self_e, double &self_de) {
  // CHRIS:: WHAT DOES THIS FUNCTION DO? WHICH LAMBDAS SHOULD WE USE??
  //         And why was it not using the member variables?
  DEBUG(14, "\t\tself energy - dipole-dipole interaction");
  const double alpha = m_A_crf_lambda_n * A_alpha + m_B_crf_lambda_n*B_alpha;

  DEBUG(15, "\t\t\talpha(A): " << A_alpha * math::four_pi_eps_i << " alpha(B): " << B_alpha * math::four_pi_eps_i << " alpha(lambda): " << alpha * math::four_pi_eps_i);

  const double d_alpha = m_n * (m_A_crf_lambda_n_1 * A_alpha - m_B_crf_lambda_n_1 * B_alpha);

  DEBUG(15, "\t\t\tn(): " << n() << " A_crf_lambda_n_1: " << A_crf_lambda_n_1() << " B_crf_lambda_n_1: " << B_crf_lambda_n_1());
  DEBUG(15, "\t\t\td_alpha(lambda)/dlambda:  " << d_alpha * math::four_pi_eps_i);

  self_e = 0.5 * alpha * e_i2;
  self_de = 0.5 * e_i2 *d_alpha;
}

/**
 * helper function to calculate the self energy 
 * at a given atom (damped).
 */
inline void interaction::Perturbed_Nonbonded_Term
::self_energy_soft_interaction(double A_alpha, double B_alpha, double e_i2,
        double A_e_0, double B_e_0, double p,
        double &self_e, double &self_de) {

  DEBUG(14, "\t\tself energy - dipole-dipole interaction");
  const double alpha = m_A_crf_lambda_n * A_alpha + m_B_crf_lambda_n * B_alpha;
  const double d_alpha = m_n * (m_A_crf_lambda_n_1 * A_alpha -
          m_B_crf_lambda_n_1 * B_alpha);
  const double e_0 = m_A_crf_lambda_n * A_e_0 + m_B_crf_lambda_n * B_e_0;
  //const double d_e_0 = m_n * (m_A_crf_lambda_n_1 * A_e_0 - m_B_crf_lambda_n_1 * B_e_0);

  const double e_02 = e_0 * e_0;
  if (e_i2 <= e_02) {
    self_e = 0.5 * alpha * e_i2;
    self_de = 0.5 * e_i2 *d_alpha;
  } else {
    const double e_i = sqrt(e_i2);
    const double e_0_div_e_i = e_0 / e_i;
    const double p_minus_1 = p - 1;
    const double p_times_p = p*p;
    //const double p_plus_1 = p + 1;
    //const double e_0_mul_d_alpha = e_0 * d_alpha;
    //const double alpha_mul_d_e_0 = alpha * d_e_0;

    self_e = 0.5 * alpha * e_02 +
            alpha * e_02 / (p * (p_minus_1)) *
            (-p_times_p +
            (e_i / e_0)*(p_times_p - 1) +
            pow(e_0_div_e_i, p_minus_1)
            );
    self_de = d_alpha *  e_02 *
            (0.5 + 1/(p * p_minus_1) * (-p_times_p + (p_times_p-1) *
            (e_i/e_0) + pow(e_0 / e_i, p_minus_1)));

  }
}



