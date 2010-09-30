/**
 * @file eds_nonbonded_term.cc
 * inline methods of Eds_Nonbonded_Term
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

/**
 * helper function to initialize the constants.
 */
inline void interaction::Eds_Nonbonded_Term
::init(simulation::Simulation const &sim) {
  switch(sim.param().force.interaction_function){
    case simulation::lj_crf_func :
      
      // Force
      m_cut3i =
              1.0 / ( sim.param().nonbonded.rf_cutoff
              * sim.param().nonbonded.rf_cutoff
              * sim.param().nonbonded.rf_cutoff);
      
      m_crf = 2*(sim.param().nonbonded.epsilon - sim.param().nonbonded.rf_epsilon) *
              (1.0 + sim.param().nonbonded.rf_kappa * sim.param().nonbonded.rf_cutoff) -
              sim.param().nonbonded.rf_epsilon * (sim.param().nonbonded.rf_kappa  *
              sim.param().nonbonded.rf_cutoff *
              sim.param().nonbonded.rf_kappa  *
              sim.param().nonbonded.rf_cutoff);
      
      m_crf /= (sim.param().nonbonded.epsilon +2* sim.param().nonbonded.rf_epsilon) *
              (1.0 + sim.param().nonbonded.rf_kappa * sim.param().nonbonded.rf_cutoff) +
              sim.param().nonbonded.rf_epsilon * (sim.param().nonbonded.rf_kappa  *
              sim.param().nonbonded.rf_cutoff *
              sim.param().nonbonded.rf_kappa  *
              sim.param().nonbonded.rf_cutoff);
      m_crf_cut3i = m_crf * m_cut3i;
      
      // Energy
      m_crf_2cut3i = m_crf_cut3i / 2.0;
      
      m_crf_cut = (1 - m_crf / 2.0)
              / sim.param().nonbonded.rf_cutoff;
      break;
      
    case simulation::pol_lj_crf_func :
    case simulation::pol_off_lj_crf_func :
    case simulation::cgrain_func :
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
::eds_lj_crf_interaction(const double dist2, const double dist2i, 
                         const double dist6i, const double disti,
			 const double &c6, const double  &c12,
			 const double &q,
			 double & force, double & e_nb)
{
  DEBUG(14, "\t\tnonbonded term");
   
  // loop over eds states
  
  const double q_eps = q * math::four_pi_eps_i;
  const double c12_dist6i = c12 * dist6i;
  e_nb = (c12_dist6i - c6) * dist6i +
          q_eps * (disti - m_crf_2cut3i * dist2 - m_crf_cut);
  
  force = (c12_dist6i + c12_dist6i - c6) * 6.0 * dist6i * dist2i +
          q_eps * (disti * dist2i + m_crf_cut3i);
  DEBUG(15, "\t\tq=" << q << " 4pie=" << math::four_pi_eps_i
          << " crf_cut2i=" << m_crf_cut3i);
  

 
  
}

/**
 * helper function to calculate the force and energy for
 * the reaction field contribution for a given pair
 */
inline void interaction::Eds_Nonbonded_Term
::eds_rf_interaction(math::Vec const &r,double q,
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


inline double interaction::Eds_Nonbonded_Term
::crf_2cut3i()const
{
  return m_crf_2cut3i;
}
