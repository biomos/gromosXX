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
              const double disti,
			        const double c6, const double c12,
			        const double q,
			        double & force, double & e_nb,
                                unsigned int eps = 0,
                                const double coulomb_scaling = 1);

    /**
     * calculate the reaction field force and energy of an atom pair.
     */
    void eds_rf_interaction(math::Vec const &r, double q,
			    math::Vec & force, double & e_rf, unsigned int eps = 0);
   

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
    
    /*
     Coarse grained dielectric permittivities*/
    std::vector<double> cgrain_eps;

    /**
     * shifting parameters
     */
    double a_RFm = 0;
    double a_RFn = 0;
    
  };
  
} // interaction

// inline methods
#include "eds_nonbonded_term.cc"

#endif
