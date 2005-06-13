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
			    double const c6, double const c12,
			    double const q,
			    double & force, double & e_lj,
			    double & e_crf);

    /**
     * calculate the reaction field force and energy of an atom pair.
     */
    void rf_interaction(math::Vec const &r, double const q,
			math::Vec & force, double & e_rf);
    /**
     * calculate the force and energy of an atom pair (coarse grain).
     */
    void cgrain_interaction(math::Vec const &r,
			    double const c6, double const c12,
			    double const q,
			    double &force, double &e_lj, double &e_crf); 
    /**
     * calculate the hessian of the lj crf term.
     */
    void lj_crf_hessian(math::Vec const &r,
			double const c6, double const c12,
			double const q,
			math::Matrix &hess);
    
    /**
     * a constant.
     */
    double const crf_2cut3i()const;

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
