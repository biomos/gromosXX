/**
 * @file eds_nonbonded_term.h
 * the eds-nonbonded terms.
 */

#ifndef INCLUDED_EDS_NONBONDED_TERM_H
#define INCLUDED_EDS_NONBONDED_TERM_H

namespace interaction
{
  /**
   * @class Nonbonded_Term
   * nonbonded terms.
   */
  class Eds_Nonbonded_Term
  {
  public:
    /**
     * Constructor.
     */
    Eds_Nonbonded_Term(){};
    
    /**
     * get the coulomb constant
     */
    double crf()const {return m_crf[0];}
    
    /**
     * initialize constants
     */
    void init(simulation::Simulation const &sim);

    /**
     * calculate the force and energy of an atom pair.
     */
    void eds_lj_crf_interaction(const double dist2, const double dist6, 
			        const double &c6, const double &c12,
			        const double &q,
                                double const alpha_lj,
			        double const alpha_crf,
			        double & force, double & e_nb,
                                unsigned int eps = 0);
    /**
     * calculate the force, energy and derivative of an atom pair
     */
    void eds_pert_lj_crf_interaction(const double dist2, const double dist6, 
				     const double &A_c6, const double &A_c12,
				     const double &B_c6, const double &B_c12,
				     const double &A_q,  const double &B_q,
				     double const alpha_lj, double const alpha_crf,
				     double & force, double & e_lj, double & e_crf, 
				     double & de_lj, double & de_crf, unsigned int eps=0);
    
    /**
     * calculate the reaction field force and energy of an atom pair.
     */
    void eds_rf_interaction(math::Vec const &r, double q, double const alpha_crf,
			    math::Vec & force, double & e_rf, unsigned int eps = 0);
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
     * reaction field cutoff.
     */
    double m_rrf;
 
    
    /**
     * Energy:
     * (1-coulomb reaction field constant / 2 and
     * divided by reaction field cutoff.
     */
    std::vector<double> m_crf_cut;
    
    /*
     Coarse grained dielectric permittivities*/
    std::vector<double> cgrain_eps;

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

  };
  
} // interaction

// inline methods
#include "eds_nonbonded_term.cc"

#endif
