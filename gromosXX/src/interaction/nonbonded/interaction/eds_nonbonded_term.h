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
    double crf()const {return m_crf;}
    
    /**
     * initialize constants
     */
    void init(simulation::Simulation const &sim);

    /**
     * calculate the force and energy of an atom pair.
     */
    void eds_lj_crf_interaction(const double dist2, const double dist2i, 
                                const double dist6i, const double disti,
			        const double &c6, const double &c12,
			        const double &q,
			        double & force, double & e_nb);

    /**
     * calculate the reaction field force and energy of an atom pair.
     */
    void eds_rf_interaction(math::Vec const &r, double q,
			    math::Vec & force, double & e_rf);
    
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
    
  };
  
} // interaction

// inline methods
#include "eds_nonbonded_term.cc"

#endif
