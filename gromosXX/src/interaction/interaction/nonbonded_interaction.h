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
     * Destructor.
     */
    virtual ~nonbonded_interaction();
    /**
     * calculate the interactions.
     */
    virtual void calculate_interactions(t_simulation &simu);
    
  protected:
    /**
     * helper class to build the pairlist.
     * @TODO parametrize that one?
     * or assume an iterator and take a reference (polymorphism)
     */
    simple_pairlist<t_simulation> m_pairlist;
    
  };
  
} // interaction

// template methods
#include "nonbonded_interaction.tcc"

#endif
