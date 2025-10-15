/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file perturbed_nonbonded_term.h
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
              const double disti,
			        const double c6, const double c12,
			        const double q,
			        double & force, double & e_nb,
              unsigned int eps = 0, const double coulomb_scaling = 1);
    /**
     * calculate the force, energy and derivative of an atom pair
     */
    void eds_pert_lj_crf_interaction(const double dist2, const double dist6, 
				     const double &A_c6, const double &A_c12,
				     const double &B_c6, const double &B_c12,
				     const double &A_q,  const double &B_q,
				     double const alpha_lj, double const alpha_crf,
				     double & force, double & e_lj, double & e_crf, 
				     double & de_lj, double & de_crf, unsigned int eps=0,
             const double coulomb_scaling = 1);
    
    /**
     * calculate the reaction field force and energy of an atom pair.
     */
    void eds_rf_interaction(math::Vec const &r, double q,
			    math::Vec & force, double & e_rf, unsigned int eps = 0);

  /**
    * calculate the reaction field force and energy of an atom pair of EDS and perturbed atom.
   */
    void eds_perturbed_rf_interaction(math::Vec const &r, double const A_q, double const B_q,
        double const alpha_crf, math::Vec & force, double &e_rf, double & de_rf, bool selfterm_correction = false, unsigned int eps = 0); 
    
    /**
     * ####### ANITA ############ 
     * precalculate the energy and dh/dl of an atom pair for other lambda points.
     */
    void lj_crf_soft_interaction_ext(math::Vec const &r,
                                 double const A_c6, double const A_c12,
                                 double const B_c6, double const B_c12,
                                 double const A_q, double const B_q,
                                 double const alpha_lj, double const alpha_crf,
                                 double &A_e_lj, double &B_e_lj,
                                 double &A_e_crf, double &B_e_crf,
                                 double &A_de_lj, double &B_de_lj,
                                 double &A_de_crf, double &B_de_crf,
                                 double const lam = 0.0, unsigned int eps = 0, 
                                 double coulomb_scaling = 1);


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
    void lj_crf_scaled_interaction(const double dist2, const double dist6,
				   double const A_c6, double const A_c12,
				   double const B_c6, double const B_c12,
				   double const A_q, double const B_q,
				   double const alpha_lj, double const alpha_crf,
				   double const A_scale, double const B_scale,
				   double & force,  double &e_lj, double & e_crf, double &de_lj, 
				   double & de_crf, unsigned int eps = 0, 
                   double coulomb_scaling = 1);

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

    //ANITA
    /**
     * calculate the reaction field force and energy
     * of a perturbed atom pair (for precalclam)
     */
    void rf_soft_interaction_ext(math::Vec const &r, double const A_q, 
                                 double const B_q, double const alpha_crf, 
                                 double & A_e_rf, double & B_e_rf,
                                 double & A_de_rf, double & B_de_rf, 
                                 double const crfs_lambda, unsigned int eps=0);
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
				 double & force, double &e_lj, double & e_crf, double &de_lj, 
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

   /**
     * calculate the force and energy of an atom pair
     * with shifted RF correction and reporting extra energies
     */
    void eds_lj_shifted_crf_corr_interaction(const double dist2, const double dist6, 
              const double disti,
              const double c6, const double c12,
              const double q,
              double &force, double &e_nb,
              double &e_extra_orig, double &e_extra_phys,
              unsigned int eps = 0,
              const double coulomb_scaling = 1);

    /**
     * calculate the reaction field force and energy of an atom pair
     * with shifted RF correction and reporting extra energies
     */
    void eds_shifted_rf_corr_interaction(math::Vec const &r, double q,
          math::Vec & force, double & e_rf, 
          double &e_extra_orig, double &e_extra_phys,
          unsigned int eps = 0);
   

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

    /**
     * Energy:
     * extra energy of distance-independent part of shifted RF
     */
    std::vector<double> m_crf_cut_extra;

    /**
     * Energy:
     * reaction field constant / twice reaction field cutoff ^ 3
     */
    std::vector<double> m_crf_2cut3i;

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


    /**
     * shifting parameters
     */
    double a_RFm = 0;
    double a_RFn = 0;
    
  };
  
} // interaction

// inline methods
#include "perturbed_nonbonded_term.cc"

#endif
