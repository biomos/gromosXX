/**
 * @file nonbonded_base.h
 * the base class for the nonbonded calculations.
 */

#ifndef INCLUDED_NONBONDED_BASE_H
#define INCLUDED_NONBONDED_BASE_H

namespace interaction
{
  /**
   * @class Nonbonded_Base
   * base class of the nonbonded calculations.
   */
  class Nonbonded_Base
  {
  public:
    /**
     * Constructor.
     */
    Nonbonded_Base();
    
    /**
     * add the lj parameters for atom type i and j.
     */
    void add_lj_parameter(size_t iac_i, size_t iac_j,
			  lj_parameter_struct lj);
    
    /**
     * get the lj parameters for atom type i and j.
     */
    lj_parameter_struct const & lj_parameter(size_t iac_i, size_t iac_j);

    /**
     * set the coulomb constant 
     */
    void coulomb_constant(double const coulomb_constant);
    
    /**
     * get the coulomb constant
     */
    double coulomb_constant()const;
    
    /**
     * set the alpha for lennard jones interactions
     */
    void alpha_lj(double const alpha_lj);
    
    /**
     * get the alpha for lennard jones interactions
     */
    double alpha_lj()const;
    
    /**
     * set the alpha for coulomd interactions
     */
    void alpha_crf(double const alpha_crf);
    /**
     * get the alpha for coulomb interactions
     */
    double alpha_crf()const;
    
    /**
     * resize the lj_parameter matrix.
     */
    void resize(size_t i);

    /**
     * initialize constants
     */
    template<typename t_simulation>
    void initialize(t_simulation const &sim);
    /**
     * calculate the force and energy of an atom pair.
     */
    void lj_crf_interaction(math::Vec const &r,
			    double const c6, double const c12,
			    double const q,
			    math::Vec & force, double & e_lj,
			    double & e_crf);
    /**
     * calculate the reaction field force and energy of an atom pair.
     */
    void rf_interaction(math::Vec const &r, double const q,
			math::Vec & force, double & e_rf);
    /**
     * calculate the force, energy and dh/dl of an atom pair.
     */
    void lj_crf_soft_interaction(math::Vec const &r,
				 double const c6, double const c12,
				 double const q, double const l,
				 double const alpha_lj,
				 double const alpha_crf,
				 math::Vec & force, double &e_lj,
				 double & e_crf, double &de_lj, 
				 double & de_crf);
    /**
     * calculate the force, energy and dh/dl of an atom pair, with scaling
     */
    void lj_crf_scaled_interaction(math::Vec const &r,
				   double const c6, double const c12,
				   double const q, double const l,
				   double const alpha_lj,
				   double const alpha_crf,
				   double const scale,
				   math::Vec & force, double &e_lj,
				   double & e_crf, double &de_lj, 
				   double & de_crf);
    /**
     * calculate the reaction field force and energy
     * of a perturbed atom pair
     */
    void rf_soft_interaction(math::Vec const &r, double const q, double const l,
			     double const alpha_crf,
			     math::Vec & force, double & e_rf, double & de_rf);
    
  protected:
    /**
     * the lj parameter.
     */
    std::vector< std::vector<lj_parameter_struct> > m_lj_parameter;

    /**
     * the coulomb constant
     */
    double m_coulomb_constant;
        
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
     * Perturbation:
     * reaction field constant / 2
     */
    double m_crf_2;
    /**
     * Perturbation:
     * reaction field cutoff ^2
     */
    double m_cut2;
    

  };
  
} // interaction

// template / inline methods
#include "nonbonded_base.tcc"

#endif
