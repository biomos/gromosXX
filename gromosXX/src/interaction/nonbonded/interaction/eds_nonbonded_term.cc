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
  double crf = 0.0, crf_cut = 0.0;
  m_crf_cut.clear();
  m_crf.clear();
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
      } else { // infinity case: rf_cutoff == 0
        crf = -1;
        DEBUG(15, "nonbonded term init: m_crf: " << crf);
        m_rrf = 0.0;
        crf_cut = (1 - crf / 2.0) / sim.param().nonbonded.rf_cutoff;
      }
      m_crf_cut.push_back(crf_cut);
      m_crf.push_back(crf);
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
        double & force, double & e_nb, unsigned int eps,
        const double coulomb_scaling) {
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
    double ecrf = q_eps * (ac_disti * coulomb_scaling - (0.5 * m_crf[eps] * ac_rrfi * dist2) - m_crf_cut[eps]);
    e_nb = elj + ecrf;

    double flj = 6.0 * dist2 * dist2 * dist6i_alj_c126 * dist6i_alj_c126 
                 * (2.0 * c12_dist6i - c6);
    double fcrf = q_eps * (ac_disti * coulomb_scaling * ac_disti * ac_disti + m_crf[eps] * ac_rrfi);
    force = flj + fcrf;
    
    DEBUG(15, "nonbonded energy = " << (elj + ecrf) << ", force = " << (flj + fcrf));
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
