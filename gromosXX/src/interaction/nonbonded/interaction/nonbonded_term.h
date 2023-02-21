/**
 * @file nonbonded_term.h
 * the nonbonded terms.
 */

#ifndef INCLUDED_NONBONDED_TERM_H
#define INCLUDED_NONBONDED_TERM_H

namespace interaction
{
  struct KSpace_Element;
  
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
     * initialize constants
     */
    void init(simulation::Simulation const &sim
            , simulation::interaction_func_enum int_func
               = simulation::default_func);

    /**
     * calculate the force and energy of an atom pair.
     */
    void lj_crf_interaction(math::Vec const &r,
			    double c6, double c12,
			    double q,
			    double & force, double & e_lj,
			    double & e_crf,
                            unsigned int eps = 0,
                            const double coulomb_scaling = 1);

    /**
     * calculate the force and energy of an atom pair.
     * New version (distance calculation in outer loop)
     */
    void lj_crf_interaction_fast(double dist2,
			    double c6, double c12,
			    double q,
			    double & force, double & e_lj,
			    double & e_crf,
                            unsigned int eps = 0);

     /**
     * calculate the force and energy of an atom pair.
     */
    void lj_interaction(math::Vec const &r,
			    double c6, double c12,
			    double & force, double & e_lj);

    /**
     * calculate the force and energy of an atom pair
     * with the shifted RF correction and report the extra energies
     */
    void lj_shifted_crf_corr_interaction(math::Vec const &r,
        double c6, double c12,
        double q,
        double &force, double &e_lj, double &e_crf,
        double &e_extra_orig, double &e_extra_phys,
        unsigned int eps = 0, const double coulomb_scaling = 1);
    
    /**
     * calculate the force and energy of an atom pair. (lattice sum)
     */
    void lj_ls_interaction(math::Vec const &r,
			double c6, double c12, double charge,
			double & force, double & e_lj, double & e_ls);
    
    /**
     * calculate the force and energy of an excluded atom pair. (lattice sum)
     */
    void ls_excluded_interaction(math::Vec const &r,
			double charge, double & force, double & e_ls);
    
    /**
     * calculate the force and energy of an atom pair. (polarisable)
     */
    void pol_lj_crf_interaction(math::Vec const &r,
                     math::Vec const &rp1,
                     math::Vec const &rp2,
                     math::Vec const &rpp,
		     double c6, double c12,
		     double qi, double qj, double cgi, double cgj,
		     double f[], double &e_lj, double &e_crf,
                     unsigned int eps = 0);

    /**
     * calculate the force and energy of an atom pair. (polarisable + off atom)
     */
    void pol_off_lj_crf_interaction(math::Vec const &r,
                     math::Vec const &rm,
                     math::Vec const &rp1,
                     math::Vec const &rp2,
                     math::Vec const &rpp,
                     double c6, double c12,
                     double qi, double qj, double cgi, double cgj,
                     double f[], double &e_lj, double &e_crf,
                     unsigned int eps = 0);

    /**
     * calculate the reaction field force and energy of an atom pair.
     */
    void rf_interaction(math::Vec const &r, double q,
			math::Vec & force, double & e_rf, unsigned int eps = 0);

    /**
     * calculate the shifted reaction field force and energy of an atom pair
     * and report the extra energies
     */
    void shifted_rf_corr_interaction(math::Vec const &r, double q,
        math::Vec &force, double &e_crf,
        double &e_extra_orig, double &e_extra_phys,
        unsigned int eps = 0);
    
    /**
     * helper function to calculate the force and energy for
     * the reaction field contribution for a given pair
     * with polarisation
     */
    inline void pol_rf_interaction(math::Vec const &r,
                 math::Vec const &rp1,
                 math::Vec const &rp2,
                 math::Vec const &rpp,
                 double qi, double qj, 
                 double cgi, double cgj,
		 math::VArray &force, double &e_crf, unsigned int eps = 0);
    /**
     * helper function to calculate the self reaction field term
     */
    inline void pol_rf_self_interaction(double const qi,double &e_crf,unsigned int eps=0);


    /**
     * calculate a term of the electric field for polarisation
     */
    void electric_field_interaction(math::Vec const &r,
			      math::Vec const &rprime,
			      double qj, double charge,
			      math::Vec &e_el, unsigned int eps = 0);

    /**
     * calculate the self energy - dipole-dipole interaction (polarisation)
     */
    void self_energy_interaction(double alpha, double e_i2, double &self_e);
    
    /**
     * calculate the damped self energy - dipole-dipole interaction (polarisation)
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
			math::Matrix &hess, unsigned int eps = 0);
    
    /**
     * a constant.
     */
    double crf_2cut3i()const;
    
    /**
     * get the coulomb constant
     */
    double crf() const;
    
    /**
     * a constant.
     */
    double crf_2cut3i(int eps) const;
    
    /**
     * a constant.
     */
    double crf_cut(int eps) const;


    /**
     * calculate the force and energy of an atom pair (sasa).
     */
    inline void sasa_interaction(math::Vec const &r, double bij,
                                 double pij, double p_i, 
                                 double surface, double & e_sasa);

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
     * reaction field constant / twice reaction field cutoff ^ 3
     */
    std::vector<double> m_crf_2cut3i;
    
    /**
     * Energy:
     * (1-coulomb reaction field constant / 2 and
     * divided by reaction field cutoff.
     */
    std::vector<double> m_crf_cut;

    /**
     * Energy:
     * extra energy of shifted RF for distance-independent part
     */
    std::vector<double> m_crf_cut_extra;
    
    /**
     * Coarse grain variables
     */
    double A_cg12, A_cg6, A_cg1;
    double B_cg12, B_cg6, B_cg1;
    double C_cg12, C_cg6, C_cg1;
    std::vector<double> cgrain_eps;

    /**
     * lattice sum variables
     */
    int charge_shape;
    double charge_width_i;

    /**
     * the cutoff squared
     */
    double m_cut2;

    /**
     * shifting parameters
     */
    double a_RFm = 0;
    double a_RFn = 0;
  };



} // interaction

// inline methods
#include "nonbonded_term.cc"

#endif
