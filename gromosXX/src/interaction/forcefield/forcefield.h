/**
 * @file forcefield.h
 * the forcefield class.
 */

#ifndef INCLUDED_FORCEFIELD_H
#define INCLUDED_FORCEFIELD_H

/**
 * @namespace interaction
 * namespace that contains the classes to
 * handle the interactions between the particles.
 * (energies, forces).
 */
namespace interaction
{
  /**
   * @class forcefield
   * contains the specific interactions.
   */
  template<typename t_simulation>
  class forcefield
  {
  public:
    /**
     * Constructor
     */
    forcefield();
    /**
     * Destructor
     */
    ~forcefield();
    /**
     * add an interaction
     */
    void add_interaction(Interaction<t_simulation> *inter);
    /**
     * calculate all interactions.
     */
    void calculate_interactions(t_simulation &simu);

  protected:
    /**
     * the interactions
     */
    std::vector<Interaction<t_simulation> *> m_interaction;
  };
  
} // interaction

// inline functions
#include "forcefield.tcc"

#endif
