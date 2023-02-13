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
  double crf = 0.0, crf_cut = 0.0, crf_cut_extra = 0.0, crf_2cut3i = 0.0;
  m_crf_cut.clear();
  m_crf.clear();
  m_crf_2cut3i.clear();
  cgrain_eps.clear();
  switch (sim.param().force.interaction_function) {
    case simulation::lj_crf_func:
    case simulation::lj_shifted_crf_corr_func:
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
        if(sim.param().nonbonded.use_shift){
          a_RFm = sim.param().nonbonded.a_RFm;
          a_RFn = sim.param().nonbonded.a_RFn;
          crf_cut_extra = sim.param().nonbonded.a_RFm * pow(sim.param().nonbonded.rf_cutoff,sim.param().nonbonded.m_crf) + sim.param().nonbonded.a_RFn * pow(sim.param().nonbonded.rf_cutoff,sim.param().nonbonded.n_crf);
          crf_cut += crf_cut_extra;
        }
      } else { // infinity case: rf_cutoff == 0
        crf = -1;
        DEBUG(15, "nonbonded eds term init: m_crf: " << crf);
        m_rrf = 0.0;

        crf_cut = (1 - crf / 2.0) / sim.param().nonbonded.rf_cutoff;
        crf_cut_extra = sim.param().nonbonded.a_RFm * pow(sim.param().nonbonded.rf_cutoff,sim.param().nonbonded.m_crf) + sim.param().nonbonded.a_RFn * pow(sim.param().nonbonded.rf_cutoff,sim.param().nonbonded.n_crf);
        crf_cut += crf_cut_extra;
      }
      crf_2cut3i = crf / (sim.param().nonbonded.rf_cutoff*sim.param().nonbonded.rf_cutoff*sim.param().nonbonded.rf_cutoff) / 2.0;
      DEBUG(15, "nonbonded eds term init: m_crf_2cut3i: " << crf_2cut3i);
      m_crf_cut.push_back(crf_cut);
      m_crf.push_back(crf);
      m_crf_cut_extra.push_back(crf_cut_extra);
      m_crf_2cut3i.push_back(crf_2cut3i);
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
        crf_2cut3i = crf / (sim.param().nonbonded.rf_cutoff*sim.param().nonbonded.rf_cutoff*sim.param().nonbonded.rf_cutoff) / 2.0;
        DEBUG(15, "nonbonded eds term init: m_crf_2cut3i: " << crf_2cut3i);
        m_crf.push_back(crf);
        m_crf_cut.push_back(crf_cut);
        m_crf_2cut3i.push_back(crf_2cut3i);
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
        const double disti,
        const double c6, const double c12,
        const double q,
        double & force, double & e_nb, unsigned int eps,
        const double coulomb_scaling) {
    DEBUG(14, "\t\tnonbonded term");

    const double dist6i = 1.0 / (dist6);
    const double q_eps = q * math::four_pi_eps_i;
    const double c12_dist6i = c12 * dist6i;
    const double dist2i = 1.0 / dist2;
    
    double elj = (c12_dist6i - c6) * dist6i;
    double ecrf = q_eps * (disti * coulomb_scaling - m_crf_2cut3i[eps] * dist2 - m_crf_cut[eps]);
    e_nb = elj + ecrf;

    double flj = (c12_dist6i + c12_dist6i - c6) * 6.0 * dist6i * dist2i;
    double fcrf = q_eps * (disti * coulomb_scaling * dist2i + 2.0 * m_crf_2cut3i[eps]);
    force = flj + fcrf;
    
    DEBUG(15, "nonbonded energy = " << (elj + ecrf) << ", force = " << (flj + fcrf));
}

/**
 * helper function to calculate the force and energy for
 * a given atom pair with the shifted RF, reporting the extra energies
 */
inline void interaction::Eds_Nonbonded_Term
::eds_lj_shifted_crf_corr_interaction(const double dist2, const double dist6, 
        const double disti,
        const double c6, const double c12,
        const double q,
        double & force, double & e_nb,
        double &e_extra_orig, double &e_extra_phys,
        unsigned int eps,
        const double coulomb_scaling) {
    DEBUG(14, "\t\tnonbonded term");

    const double dist4 = dist2 * dist2;
    const double dist6i = 1.0 / (dist6);
    const double q_eps = q * math::four_pi_eps_i;
    const double c12_dist6i = c12 * dist6i;
    const double dist2i = 1.0 / dist2;
    
    double elj = (c12_dist6i - c6) * dist6i;
    double ecrf = q_eps * (disti * coulomb_scaling - m_crf_2cut3i[eps] * dist2 + a_RFm * dist4 + a_RFn * dist6 - m_crf_cut[eps]);
    e_extra_orig = - q_eps * (a_RFm * dist4 + a_RFn * dist6 - m_crf_cut_extra[eps]);
    e_extra_phys = - q_eps * (a_RFm * dist4 + a_RFn * dist6 - m_crf_cut[eps]);
    e_nb = elj + ecrf;

    double flj = (c12_dist6i + c12_dist6i - c6) * 6.0 * dist6i * dist2i;
    double fcrf = q_eps * (disti * coulomb_scaling * dist2i - 4 * a_RFm * dist2 - 6 * a_RFn * dist4 + 2.0 * m_crf_2cut3i[eps]);
    force = flj + fcrf;
    
    DEBUG(15, "nonbonded energy = " << (elj + ecrf) << ", force = " << (flj + fcrf));
}

/**
 * helper function to calculate the force and energy for
 * the reaction field contribution for a given pair
 */
inline void interaction::Eds_Nonbonded_Term
::eds_rf_interaction(math::Vec const &r, double q, 
		     math::Vec &force, double &e_crf, unsigned int eps)
{
  const double dist2 = abs2(r);
  
  force = q * math::four_pi_eps_i * 2.0 * m_crf_2cut3i[eps] * r;

  e_crf = q * math::four_pi_eps_i * ( -m_crf_2cut3i[eps] * dist2 - m_crf_cut[eps]);
  DEBUG(15, "dist2 " << dist2 );
  DEBUG(15, "q*q   " << q );
  DEBUG(15, "rf energy = " << e_crf);
  
}

/**
 * helper function to calculate the force and energy for
 * the reaction field contribution for a given pair
 * with shifted RF and corrected energies
 */
inline void interaction::Eds_Nonbonded_Term
::eds_shifted_rf_corr_interaction(math::Vec const &r, double q, 
         math::Vec &force, double &e_crf, 
         double &e_extra_orig, double &e_extra_phys,
         unsigned int eps)
{
  const double dist2 = abs2(r);
  const double dist4 = dist2 * dist2;
  const double dist6 = dist2 * dist4;
  
  force = q * math::four_pi_eps_i * (-4 * a_RFm * dist2 - 6 * a_RFn * dist4 + 2.0 * m_crf_2cut3i[eps]) * r;

  e_crf = q * math::four_pi_eps_i * ( -m_crf_2cut3i[eps] * dist2 + a_RFm * dist4 + a_RFn * dist6 - m_crf_cut[eps]);
  e_extra_orig = - q * math::four_pi_eps_i * (a_RFm * dist4 + a_RFn * dist6 - m_crf_cut_extra[eps]);
  e_extra_phys = - q * math::four_pi_eps_i * (a_RFm * dist4 + a_RFn * dist6 - m_crf_cut[eps]);

  DEBUG(15, "dist2 " << dist2 );
  DEBUG(15, "q*q   " << q );
  DEBUG(15, "rf energy = " << e_crf);
  
}