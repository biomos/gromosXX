/**
 * @file perturbed_nonbonded_term.h
 * perturbed nonbonded terms.
 */

#ifndef INCLUDED_PERTURBED_NONBONDED_TERM_H
#define INCLUDED_PERTURBED_NONBONDED_TERM_H

namespace interaction
{
  /**
   * @class Perturbed_Nonbonded_Term
   * base class of the nonbonded calculations.
   */
  class Perturbed_Nonbonded_Term
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Nonbonded_Term(){};
    
    /**
     * get the coulomb constant
     */
    double crf()const {return m_crf[0];}
    
    /**
     * initialize constants
     */
    void init(simulation::Simulation const &sim);

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
				 double & de_crf,
                                 unsigned int eps = 0);
    
    /**
     * calculate the force, energy and dh/dl of an atom pair. for polarisation
     */
    void pol_lj_crf_soft_interaction(math::Vec const &r, math::Vec const &rp1,
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
                                     double &de_lj, double &de_crf,
                                     unsigned int eps = 0);
    /**
     * calculate the force, energy and dh/dl of an atom pair. for polarisation + off atom
     */
    void pol_off_lj_crf_soft_interaction(math::Vec const &r, math::Vec const &rm, math::Vec const &rp1,
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
                                     double &de_lj, double &de_crf,
                                     unsigned int eps = 0);


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
				   double & de_crf,
                                   unsigned int eps = 0);

    /**
     * calculate the reaction field force and energy
     * of a perturbed atom pair
     */
    void rf_soft_interaction(math::Vec const &r, 
			     double const A_q, double const B_q, 
			     double const alpha_crf,
			     math::Vec & force, double & e_rf,
			     double & de_rf,
			     bool selfterm_correction = false,
                             unsigned int eps = 0);
    /**
     * calculate the reaction field force and energy
     * of a perturbed atom pair (with polarisation)
     */
    void pol_rf_soft_interaction(math::Vec const &r,
                             math::Vec const &rp1,
                             math::Vec const &rp2,
                             math::Vec const &rpp,
			     double const A_qi, double const A_qj,
                             double const B_qi, double const B_qj,
                             double cqi, double cqj,
			     double const alpha_crf,
			     double force[], double & e_rf,
			     double & de_rf,
			     bool selfterm_correction = false,
                             unsigned int eps = 0);

    /**
     * calculate the self reaction field energy with polarisation
     */
    void pol_rf_self_soft_interaction(double const A_qi,double const B_qi,
	                              double & e_rf, double & de_rf,
                                      bool selfterm_correction = false,
                                      unsigned int eps = 0); 

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
			       math::Vec & force, double & e_rf,
			       double & de_rf,
			       bool selfterm_correction = false,
                               unsigned int eps = 0);

    /**
     * calculate the force, energy and dh/dl of an atom pair (coarse grain)
     */
    void cgrain_soft_interaction(math::Vec const &r,
				 double const A_c6, double const A_c12,
				 double const B_c6, double const B_c12,
				 double const A_q, double const B_q,
				 double const alpha_lj,
				 double const alpha_crf,
				 double & force1, double & force6, double & force12,
				 double &e_lj, double & e_crf, double &de_lj, 
				 double & de_crf);
    
    /**
     * calculate the perturbed electric field term.
     */
    void electric_field_soft_interaction(math::Vec const &r, 
                       math::Vec const &rprime, 
		       double const alpha_crf,
                       double A_qj, double B_qj, double cgj, 
		       math::Vec &e_el,
                       unsigned int eps = 0);
    
    /**
     * calculate the self energy - dipole-dipole interaction (polarisation)
     */
    void self_energy_soft_interaction(double A_alpha, double B_alpha, 
                                      double e_i2, double &self_e, double &self_de);
    
    /**
     * calculate the damped self energy - dipole-dipole interaction (polarisation)
     */
    void self_energy_soft_interaction(double A_alpha, double B_alpha, double e_i2, 
                                      double A_e_0, double B_e_0, double p,
                                      double &self_e, double &self_de);

    /**
     * Perturbation:
     * lambda value for lj interaction for state A
     */
    double const & A_lj_lambda()const;
    /**
     * Perturbation:
     * lambda value for crf interaction for state A
     */
    double const & A_crf_lambda()const;
    /**
     * Perturbation:
     * lambda value for lj interaction for state B
     */
    double const & B_lj_lambda()const;
    /**
     * Perturbation:
     * lambda value for crf interaction for state B
     */
    double const & B_crf_lambda()const;
    /**
     * Perturbation:
     * lambda exponent:
     */
    int const & n()const;
    /**
     * Perturbation:
     * lambda value for lj interaction for state A to the power nlam
     */
    double const & A_lj_lambda_n()const;
    /**
     * Perturbation:
     * lambda value for crf interaction for state A to the power nlam
     */
    double const & A_crf_lambda_n()const;
    /**
     * Perturbation:
     * lambda value for lj interaction for state B to the power nlam
     */
    double const & B_lj_lambda_n()const;
    /**
     * Perturbation:
     * lambda value for crf interaction for state B to the power nlam
     */
    double const & B_crf_lambda_n()const;
    /**
     * Perturbation:
     * lambda value for lj interaction for state A to the power nlam-1
     */
    double const & A_lj_lambda_n_1()const;
    /**
     * Perturbation:
     * lambda value for crf interaction for state A to the power nlam-1
     */
    double const & A_crf_lambda_n_1()const;   
    /**
     * Perturbation:
     * lambda value for lj interaction for state B to the power nlam-1
     */
    double const & B_lj_lambda_n_1()const;
     /**
     * Perturbation:
     * lambda value for crf interaction for state B to the power nlam-1
     */
    double const & B_crf_lambda_n_1()const;
    /**
     * Perturbation:
     * set the lambdas
     */
    void set_lambda(double const lj_lambda, double const ljs_lambda,
                    double const crf_lambda, double const crfs_lambda,
                    double const lj_lambda_derivative, double const ljs_lambda_derivative,
                    double const crf_lambda_derivative, double const crfs_lambda_derivative,
                    int const n);
    
  protected:
    /**
     * Force:
     * reaction field constant.
     */
    std::vector<double> m_crf;

   /**
     * Force:
     * inverse reaction field cutoff to the power of 3.
     */
    std::vector<double> m_cut3i;

    /**
     * Force:
     * coulomb reaction field constant divided
     * by the reaction field cutoff to the power of 3.
     */
    std::vector<double> m_crf_cut3i;

    /**
     * Energy:
     * (1-coulomb reaction field constant / 2 and
     * divided by reaction field cutoff.
     */
    std::vector<double> m_crf_cut;
    /**
     * Energy:
     * reaction field constant / 2
     */
    std::vector<double> m_crf_2;
    /**
     * Perturbation:
     * reaction field cutoff ^2
     */
    double m_cut2;
    /**
     * Perturbation:
     * lambda value for lj interaction for state A
     */
    double m_A_lj_lambda;
     /**
     * Perturbation:
     * lambda value for lj softness for state A
     */
    double m_A_ljs_lambda;
    /**
     * Perturbation:
     * lambda value for crf interaction for state A
     */
    double m_A_crf_lambda;
     /**
     * Perturbation:
     * lambda value for crf softness for state A
     */
    double m_A_crfs_lambda;
    /**
     * Perturbation
     * lambda value for lj softness for state A squared.
     */
    double m_A_ljs_lambda2;
    /**
     * Perturbation
     * lambda value for crf softness for state A squared.
     */
    double m_A_crfs_lambda2;
    /**
     * Perturbation
     * lambda value for lj interaction for state A to the power nlam
     */
    double m_A_lj_lambda_n;
    /**
     * Perturbation
     * lambda value for crf interaction for state A to the power nlam
     */
    double m_A_crf_lambda_n;
    /**
     * Perturbation
     * lambda value for lj interaction for state A to the power nlam-1
     */
    double m_A_lj_lambda_n_1;
    /**
     * Perturbation
     * lambda value for crf interaction for state A to the power nlam-1
     */
    double m_A_crf_lambda_n_1;
    /**
     * Perturbation:
     * lambda value for lj interaction for state B
     */
    double m_B_lj_lambda;
     /**
     * Perturbation:
     * lambda value for lj softness for state B
     */
    double m_B_ljs_lambda;
    /**
     * Perturbation:
     * lambda value for crf interaction for state B
     */
    double m_B_crf_lambda;
     /**
     * Perturbation:
     * lambda value for crf softness for state B
     */
    double m_B_crfs_lambda;
    /**
     * Perturbation
     * lambda value for lj softness for state B squared.
     */
    double m_B_ljs_lambda2;
    /**
     * Perturbation
     * lambda value for crf softness for state B squared.
     */
    double m_B_crfs_lambda2;
    /**
     * Perturbation
     * lambda value for lj interaction for state B to the power nlam
     */
    double m_B_lj_lambda_n;
    /**
     * Perturbation
     * lambda value for crf interaction for state B to the power nlam
     */
    double m_B_crf_lambda_n;
    /**
     * Perturbation
     * lambda value for lj interaction for state B to the power nlam-1
     */
    double m_B_lj_lambda_n_1;
    /**
     * Perturbation
     * lambda value for crf interaction for state B to the power nlam-1
     */
    double m_B_crf_lambda_n_1;
    /**
     * Perturbation:
     * exponent to lambda.
     */
    double m_lambda_exp;
    /**
     * Perturbation
     * exponent to lambda as an int
     */
    int m_n;
    /**
     * Coarse grain variables
     */
    double A_cg12, A_cg6, A_cg1;
    double B_cg12, B_cg6, B_cg1;
    double A_C_cg12, A_C_cg6, A_C_cg1;
    double B_C_cg12, B_C_cg6, B_C_cg1;
    std::vector<double> cgrain_eps;
    double rs_lj, rs_c;
    double nb_cutoff;
  };
  
} // interaction

// template / inline methods
#include "perturbed_nonbonded_term.cc"

#endif
