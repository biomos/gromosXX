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
   * @class nonbonded_interaction
   * calculates the nonbonded interactions.
   */
  template<typename t_simulation>
  class nonbonded_interaction : public interaction<t_simulation>
  {
  public:
    /**
     * @struct lj_parameter_struct
     * Lennard Jones interaction parameter.
     */
    struct lj_parameter_struct
    {
      double c6;
      double c12;
      double cs6;
      double cs12;
    };
    
    /**
     * Destructor.
     */
    virtual ~nonbonded_interaction();

    /**
     * calculate the interactions.
     */
    virtual void calculate_interactions(t_simulation &simu);

    /**
     * add the lj parameters for atom type i and j.
     */
    void add_lj_parameter(int i, int j, lj_parameter_struct lj);
    
    /**
     * resize the lj_parameter matrix.
     */
    void resize(size_t i);
    
  protected:
    /**
     * helper class to build the pairlist.
     * @TODO parametrize that one?
     * or assume an iterator and take a reference (polymorphism)
     */
    simple_pairlist<t_simulation> m_pairlist;

    std::vector< std::vector<lj_parameter_struct> > m_lj_parameter;

  };
  
} // interaction

// template methods
#include "nonbonded_interaction.tcc"

#endif
