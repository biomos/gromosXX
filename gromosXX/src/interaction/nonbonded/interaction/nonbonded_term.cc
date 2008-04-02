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
  switch(sim.param().force.interaction_function){
  case simulation::lj_crf_func :
  case simulation::pol_lj_crf_func :
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
    break;
    
  case simulation::cgrain_func :
    // cgrain
    A_cg12= - (12.0 * (12 + 4)) / (pow(sim.param().longrange.rf_cutoff, 12 + 3));
    A_cg6=  - (6.0  * (6  + 4)) / (pow(sim.param().longrange.rf_cutoff, 6  + 3));
    A_cg1=  - (1.0  * (1  + 4)) / (pow(sim.param().longrange.rf_cutoff, 1  + 3));     
    
    B_cg12=   (12.0 * (12 + 3)) / (pow(sim.param().longrange.rf_cutoff, 12 + 4));
    B_cg6=    (6.0  * (6  + 3)) / (pow(sim.param().longrange.rf_cutoff, 6  + 4));
    B_cg1=    (1.0  * (1  + 3)) / (pow(sim.param().longrange.rf_cutoff, 1  + 4));     

    C_cg12=   ((12 + 3) * (12 + 4)) / (12.0 * pow(sim.param().longrange.rf_cutoff, 12));
    C_cg6=    ((6  + 3) * (6  + 4)) / (12.0 * pow(sim.param().longrange.rf_cutoff, 6 ));
    C_cg1=    ((1  + 3) * (1  + 4)) / (12.0 * pow(sim.param().longrange.rf_cutoff, 1 ));
    
    cgrain_eps = sim.param().cgrain.EPS;
    break;
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
inline void interaction::Nonbonded_Term
::lj_crf_interaction(math::Vec const &r,
		     double c6, double c12,
		     double q,
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
 * a given atom pair (polarizable).
 */
inline void interaction::Nonbonded_Term
::pol_lj_crf_interaction(math::Vec const &r,
                     math::Vec const &rp1,
                     math::Vec const &rp2,
                     math::Vec const &rpp,
		     double c6, double c12,
		     double qi, double qj, double cgi, double cgj,
		     std::vector<double> &f, double &e_lj, double &e_crf)
{
  DEBUG(14, "\t\tnonbonded term");
  
  assert(abs2(r) != 0);

  const double dist2 = abs2(r);
  const double dist2p1 = abs2(rp1);
  const double dist2p2 = abs2(rp2);
  const double dist2pp = abs2(rpp);
  const double dist2i = 1.0 / dist2;
  const double dist2p1i = 1.0 / dist2p1;
  const double dist2p2i = 1.0 / dist2p2;
  const double dist2ppi = 1.0 / dist2pp;

  const double dist6i = dist2i * dist2i * dist2i;

  const double disti = sqrt(dist2i);
  const double distp1i = sqrt(dist2p1i);
  const double distp2i = sqrt(dist2p2i);
  const double distppi = sqrt(dist2ppi);

  const double c12_dist6i = c12 * dist6i;
  const double eps = math::four_pi_eps_i;
  const double q_eps = (qi-cgi)*(qj-cgj) * eps;
  const double q_epsp1 = (qi-cgi)*cgj * eps;
  const double q_epsp2 = cgi*(qj-cgj) * eps;
  const double q_epspp = cgi*cgj * eps;
  
  e_lj = (c12_dist6i - c6) * dist6i;

  e_crf = q_eps * (disti - m_crf_2cut3i * dist2 - m_crf_cut)
    + q_epsp1 * (distp1i - m_crf_2cut3i * dist2p1 - m_crf_cut)
    + q_epsp2 * (distp2i - m_crf_2cut3i * dist2p2 - m_crf_cut)
    + q_epspp * (distppi - m_crf_2cut3i * dist2pp - m_crf_cut);

  f[0] = (c12_dist6i + c12_dist6i - c6) * 6.0 * dist6i * dist2i + 
         q_eps * (disti * dist2i + m_crf_cut3i);
  f[1] = q_epsp1 * (distp1i * dist2p1i + m_crf_cut3i);
  f[2] = q_epsp2 * (distp2i * dist2p2i + m_crf_cut3i);
  f[3] = q_epspp * (distppi * dist2ppi + m_crf_cut3i);

  DEBUG(15, "\t\tq=" << qi*qj << " 4pie=" << math::four_pi_eps_i 
	<< " crf_cut2i=" << m_crf_cut3i);
  
}

/**
 * helper function to calculate the force and energy for
 * the reaction field contribution for a given pair
 */
inline void interaction::Nonbonded_Term
::rf_interaction(math::Vec const &r,double q,
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

/**
 * helper function to calculate the force and energy for
 * the reaction field contribution for a given pair
 * with polarization
 */
inline void interaction::Nonbonded_Term
::pol_rf_interaction(math::Vec const &r,
                 math::Vec const &rp1,
                 math::Vec const &rp2,
                 math::Vec const &rpp,
                 double qi, double qj, 
                 double cgi, double cgj,
		 math::VArray &force, double &e_crf)
{
  const double dist2 = abs2(r);
  const double dist2p1 = abs2(rp1);
  const double dist2p2 = abs2(rp2);
  const double dist2pp = abs2(rpp);

  const double eps = math::four_pi_eps_i;
  const double qeps = (qi-cgi)*(qj-cgj) * eps;
  const double qepsp1 = (qi-cgi)*cgj * eps;
  const double qepsp2 = cgi*(qj-cgj) * eps;
  const double qepspp = cgi*cgj * eps;
  
  force(0) = qeps *  m_crf_cut3i * r;
  force(1) = qepsp1 *  m_crf_cut3i * rp1;
  force(2) = qepsp2 *  m_crf_cut3i * rp2;
  force(3) = qepspp *  m_crf_cut3i * rpp;

  e_crf = qeps*( -m_crf_2cut3i*dist2 - m_crf_cut)
          + qepsp1*( -m_crf_2cut3i*dist2p1 - m_crf_cut)
          + qepsp2*( -m_crf_2cut3i*dist2p2 - m_crf_cut)
          + qepspp*( -m_crf_2cut3i*dist2pp - m_crf_cut);
}


/**
 * helper function to calculate the force and energy for
 * a given atom pair in the coarse grain model
 */
inline void interaction::Nonbonded_Term
::cgrain_interaction(math::Vec const &r,
                     double c6, double c12,
                     double q,
                     double &force, double &e_lj, double &e_crf)
{
  assert(abs2(r) != 0);
  const double dist2 = abs2(r);
  const double dist2i = 1.0 / dist2;
  const double disti = sqrt(dist2i);
  const double dist6i = dist2i * dist2i * dist2i;
  const double dist12i = dist6i * dist6i;
  const double dist = 1.0 / disti;

  const double q_eps = q * math::four_pi_eps_i / cgrain_eps;

  // const double c12_dist6i = c12 * dist6i;

  e_crf = (q_eps * (disti
                    - A_cg1  / 3 * dist2 * dist
                    - B_cg1  / 4 * dist2 * dist2
                    - C_cg1 ));
  
  e_lj =  (c12 * (dist12i
		  - A_cg12 / 3 * dist2 * dist
		  - B_cg12 / 4 * dist2 * dist2
		  - C_cg12))
    -     (c6 *  (dist6i
		  - A_cg6  / 3 * dist2 * dist
		  - B_cg6  / 4 * dist2 * dist2 
		  - C_cg6 ));
  
  force = c12 * (12.0 * dist12i * disti + A_cg12 * dist2 + B_cg12 * dist2 * dist) * disti -
           c6 * ( 6.0 * dist6i *  disti + A_cg6  * dist2 + B_cg6  * dist2 * dist) * disti +
               q_eps * (         dist2i + A_cg1  * dist2 + B_cg1  * dist2 * dist) * disti;
  
  
  std::cout.precision(10);
  
  DEBUG(11, "r_ij= " << dist 
        << " e_lj=" << e_lj << " e_crf=" << e_crf 
        << " force=" << force);
  
}

inline double interaction::Nonbonded_Term
::crf_2cut3i()const
{
  return m_crf_2cut3i;
}

/**
 * helper function to calculate a term of the electric field 
 * at a given position for the polarization
 */
inline void interaction::Nonbonded_Term
::electric_field_interaction(math::Vec const &r, 
                       math::Vec const &rprime, 
                       double qj, double charge, 
                       math::Vec &e_el) {

  DEBUG(14, "\t\tenergy field term for polarization");

  assert(abs2(r) != 0);
  assert(abs2(rprime) != 0);
  const double distj = abs2(r);
  const double distp = abs2(rprime);
  const double distji = 1/(distj*sqrt(distj));
  const double distpi = 1/(distp*sqrt(distp));
  const double q_eps = (qj-charge) * math::four_pi_eps_i;
  const double q_epsp = charge * math::four_pi_eps_i;

  e_el = q_eps*(distji + m_crf_cut3i)*r + q_epsp*(distpi + m_crf_cut3i)*rprime;
}

/**
 * helper function to calculate the self energy 
 * at a given atom.
 */
inline void interaction::Nonbonded_Term
::self_energy_interaction(double alpha, double e_i2, double &self_e) {

  DEBUG(14, "\t\tself energy - dipole-dipole interaction");
  self_e = 0.5 * alpha * e_i2;
}

/**
 * helper function to calculate the self energy 
 * at a given atom (damped).
 */
inline void interaction::Nonbonded_Term
::self_energy_interaction(double alpha, double e_i2, double e_0, double p,
                       double &self_e) {

  DEBUG(14, "\t\tself energy - dipole-dipole interaction");
  const double e_02 = e_0 * e_0; 
  if (e_i2 <= e_02) {
    self_e = 0.5 * alpha * e_i2;
  } else {
    self_e = 0.5 * alpha * e_02 / (p - 1) *
             (p + 1 - 2 * pow(e_0 / sqrt(e_i2), p - 1));
  } 
}

inline void
interaction::Nonbonded_Term::lj_crf_hessian(math::Vec const &r,
				    double c6, double c12,
				    double q,
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

