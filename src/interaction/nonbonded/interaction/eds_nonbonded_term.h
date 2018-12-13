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
     * calculate the reaction field force and energy of an atom pair.
     */
    void eds_rf_interaction(math::Vec const &r, double q, double const alpha_crf,
			    math::Vec & force, double & e_rf, unsigned int eps = 0);
   

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
    
  };
  
} // interaction

// inline methods
#include "eds_nonbonded_term.cc"

#endif
