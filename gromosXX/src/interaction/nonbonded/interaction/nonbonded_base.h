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
    Nonbonded_Base(){};
    
    /**
     * resize the lj_parameter matrix.
     */
    void resize(size_t i);
    /**
     * the lj parameter.
     */
    std::vector<std::vector<lj_parameter_struct> > & lj_parameter()
    {
      return m_lj_parameter;
    }
    /**
     * get the lj parameters for atom type i and j.
     */
    lj_parameter_struct const & lj_parameter(size_t iac_i, size_t iac_j){
      assert(iac_i < m_lj_parameter.size());
      assert(iac_j < m_lj_parameter[iac_i].size());
      
      return m_lj_parameter[iac_i][iac_j];
    }
    
    /**
     * add the lj parameters for atom type i and j.
     */
    void add_lj_parameter(size_t iac_i, size_t iac_j,
			  lj_parameter_struct lj);
    
    /**
     * get the coulomb constant
     */
    double crf()const {return m_crf;}
    
    /**
     * initialize constants
     */
    void initialize(simulation::Simulation const &sim);
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
				 double const A_c6, double const A_c12,
				 double const B_c6, double const B_c12,
				 double const A_q, double const B_q,
				 double const alpha_lj,
				 double const alpha_crf,
				 double & force1, double & force6, double & force12,
				 double &e_lj, double & e_crf, double &de_lj, 
				 double & de_crf);

    /**
     * calculate the force, energy and dh/dl of an atom pair for
     * which the interaction is scaled.
     */
    void lj_crf_scaled_interaction(math::Vec const &r,
				   double const A_c6, double const A_c12,
				   double const B_c6, double const B_c12,
				   double const A_q, double const B_q,
				   double const alpha_lj, double const alpha_crf,
				   double const A_scale, double const B_scale,
				   double & force1, double & force6, double & force12,
				   double &e_lj, double & e_crf, double &de_lj, 
				   double & de_crf);


    /**
     * calculate the reaction field force and energy
     * of a perturbed atom pair
     */
    void rf_soft_interaction(math::Vec const &r, 
			     double const A_q, double const B_q, 
			     double const l,
			     double const alpha_crf,
			     math::Vec & force, double & e_rf, double & de_rf);

    /**
     * calculate the reaction field force and energy
     * of a perturbed atom pair for which the interaction
     * is scaled.
     */
    void rf_scaled_interaction(math::Vec const &r, 
			       double const A_q, double const B_q, 
			       double const l,
			       double const alpha_crf,
			       double const A_scale, double const B_scale,
			       math::Vec & force, double & e_rf, double & de_rf);
    
    void lj_crf_hessian(math::Vec const &r,
			double const c6, double const c12,
			double const q,
			math::Matrix &hess);

    double const crf_2cut3i()const;
     /**
     * Perturbation:
     * lambda value for state A
     */
    double const A_lambda()const;
    /**
     * Perturbation:
     * lambda value for state B
     */
    double const B_lambda()const;
    /**
     * Perturbation:
     * lambda value for state A to the power nlam
     */
    double const A_lambda_n()const;
    /**
     * Perturbation:
     * lambda value for state B to the power nlam
     */
    double const B_lambda_n()const;
    /**
     * Perturbation:
     * lambda value for state A to the power nlam-1
     */
    double const A_lambda_n_1()const;
    /**
     * Perturbation:
     * lambda value for state B to the power nlam-1
     */
    double const B_lambda_n_1()const;
    /**
     * Perturbation:
     * set the lambdas
     */
    void set_lambda(double const l, int const n);
    

  protected:
    /**
     * the lj parameter.
     */
    std::vector< std::vector<lj_parameter_struct> > m_lj_parameter;

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
     * Perturbation:
     * reaction field constant / 2
     */
    double m_crf_2;
    /**
     * Perturbation:
     * reaction field cutoff ^2
     */
    double m_cut2;
    /**
     * Perturbation:
     * lambda value for state A
     */
    double m_A_lambda;
    /**
     * Perturbation
     * A lambda squared.
     */
    double m_A_lambda2;
    /**
     * Perturbation:
     * lambda value for state B
     */
    double m_B_lambda;
    /**
     * Perturbation:
     * square lambda value for state B
     */
    double m_B_lambda2;
    /**
     * Perturbation:
     * lambda value for state A to the power nlam
     */
    double m_A_lambda_n;
    /**
     * Perturbation:
     * lambda value for state B to the power nlam
     */
    double m_B_lambda_n;
    /**
     * Perturbation:
     * lambda value for state A to the power nlam-1
     */
    double m_A_lambda_n_1;
    /**
     * Perturbation:
     * lambda value for state B to the power nlam-1
     */
    double m_B_lambda_n_1;

    /**
     * Perturbation:
     * exponent to lambda.
     */
    double m_lambda_exp;

  };
  
} // interaction

// template / inline methods
#include "nonbonded_base.tcc"

#endif
