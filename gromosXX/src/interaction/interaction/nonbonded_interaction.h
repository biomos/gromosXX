/**
 * @file nonbonded_interaction.h
 * the non bonded interactions:
 * Lennard-Jones and Coulomb interactions.
 */

#ifndef INCLUDED_NONBONDED_INTERACTION_H
#define INCLUDED_NONBONDED_INTERACTION_H

namespace interaction
{
  /**
   * @class Nonbonded_Interaction
   * calculates the nonbonded interactions.
   */
  template<typename t_simulation, typename t_pairlist>
  class Nonbonded_Interaction : public Interaction<t_simulation>
  {
  public:    
    /**
     * Constructor.
     */
    Nonbonded_Interaction();
    
    /**
     * Destructor.
     */
    virtual ~Nonbonded_Interaction();

    /**
     * calculate the interactions.
     */
    virtual void calculate_interactions(t_simulation &sim);

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
     * resize the lj_parameter matrix.
     */
    void resize(size_t i);

    /**
     * pairlist accessor
     */
    t_pairlist & pairlist();
    
  protected:
    /**
     * @enum nonbonded_type_enum
     * type of nonbonded interaction to calculate.
     */
    enum nonbonded_type_enum { shortrange, longrange };

    /**
     * helper class to build the pairlist.
     * @TODO parametrize that one?
     * or assume an iterator and take a reference (polymorphism)
     */
    t_pairlist m_pairlist;

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
     * the long-range force.
     */
    math::VArray m_longrange_force;
    
    /**
     * calculate the interactions.
     */
    virtual void do_interactions(t_simulation &sim,
				 typename t_pairlist::iterator it, 
				 typename t_pairlist::iterator to,
				 nonbonded_type_enum range = shortrange);

    /**
     * calculate the 1,4-interactions.
     */
    virtual void do_14_interactions(t_simulation &sim);

    /**
     * calculate the RF contributions for excluded atoms.
     */
    virtual void do_RF_excluded_interactions(t_simulation &sim);
    
    /**
     * initialize constants
     */
    void initialize(t_simulation const &sim);

    /**
     * calculate the force and energy of an atom pair.
     */
    void lj_crf_interaction(math::Vec const &r,
			    double const c6, double const c12,
			    double const q,
			    math::Vec & force, double & e_lj, double & e_crf);
    /**
     * calculate the reaction field force and energy of an atom pair.
     */
    void rf_interaction(math::Vec const &r, double const q,
			math::Vec & force, double & e_rf);

  };
  
} // interaction

// template methods
#include "nonbonded_interaction.tcc"

#endif
