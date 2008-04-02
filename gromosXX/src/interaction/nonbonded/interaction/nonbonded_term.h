/**
 * @file nonbonded_term.h
 * the nonbonded terms.
 */

#ifndef INCLUDED_NONBONDED_TERM_H
#define INCLUDED_NONBONDED_TERM_H

namespace interaction
{
  /**
   * @class Nonbonded_Term
   * nonbonded terms.
   */
  class Nonbonded_Term
  {
  public:
    /**
     * Constructor.
     */
    Nonbonded_Term(){};
    
    /**
     * get the coulomb constant
     */
    double crf()const {return m_crf;}
    
    /**
     * initialize constants
     */
    void init(simulation::Simulation const &sim);

    /**
     * calculate the force and energy of an atom pair.
     */
    void lj_crf_interaction(math::Vec const &r,
			    double c6, double c12,
			    double q,
			    double & force, double & e_lj,
			    double & e_crf);
    
    /**
     * calculate the force and energy of an atom pair. (polarizable)
     */
    void pol_lj_crf_interaction(math::Vec const &r,
                     math::Vec const &rp1,
                     math::Vec const &rp2,
                     math::Vec const &rpp,
		     double c6, double c12,
		     double qi, double qj, double cgi, double cgj,
		     std::vector<double> &f, double &e_lj, double &e_crf);

    /**
     * calculate the reaction field force and energy of an atom pair.
     */
    void rf_interaction(math::Vec const &r, double q,
			math::Vec & force, double & e_rf);
    
    /**
     * helper function to calculate the force and energy for
     * the reaction field contribution for a given pair
     * with polarization
     */
    inline void pol_rf_interaction(math::Vec const &r,
                 math::Vec const &rp1,
                 math::Vec const &rp2,
                 math::Vec const &rpp,
                 double qi, double qj, 
                 double cgi, double cgj,
		 math::VArray &force, double &e_crf);
    
    /**
     * calculate a term of the electric field for polarization
     */
    void electric_field_interaction(math::Vec const &r,
			      math::Vec const &rprime,
			      double qj, double charge,
			      math::Vec  &e_el);


    /**
     * calculate the self energy - dipole-dipole interaction (polarization)
     */
    void self_energy_interaction(double alpha, double e_i2, double &self_e);
    
    /**
     * calculate the damped self energy - dipole-dipole interaction (polarization)
     */
    void self_energy_interaction(double alpha, double e_i2, double e_0, double p,
                       double &self_e);

    /**
     * calculate the force and energy of an atom pair (coarse grain).
     */
    void cgrain_interaction(math::Vec const &r,
			    double c6, double c12,
			    double q,
			    double &force, double &e_lj, double &e_crf); 

    /**
     * calculate the hessian of the lj crf term.
     */
    void lj_crf_hessian(math::Vec const &r,
			double c6, double c12,
			double  q,
			math::Matrix &hess);
    
    /**
     * a constant.
     */
    double crf_2cut3i()const;

  protected:
    /**
     * Force:
     * reaction field constant.
     */
    double m_crf;
 
   /**
     * Force:
     * inverse reaction field cutoff to the power of 3.
     */
    double m_cut3i;
     
    /**
     * Force:
     * coulomb reaction field constant divided
     * by the reaction field cutoff to the power of 3.
     */
    double m_crf_cut3i;
    
    /**
     * Energy:
     * reaction field constant / twice reaction field cutoff ^ 3
     */
    double m_crf_2cut3i;
    
    /**
     * Energy:
     * (1-coulomb reaction field constant / 2 and
     * divided by reaction field cutoff.
     */
    double m_crf_cut;
    
    /**
     * Coarse grain variables
     */
    double A_cg12, A_cg6, A_cg1;
    double B_cg12, B_cg6, B_cg1;
    double C_cg12, C_cg6, C_cg1;
    double cgrain_eps;
  };
  
} // interaction

// inline methods
#include "nonbonded_term.cc"

#endif
