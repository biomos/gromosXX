/**
 * @file forcefield.h
 * the Forcefield class.
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
   * @class Forcefield
   * contains the specific interactions.
   * @TODO are the destructors called? properly?
   * clear does not call them (i guess).
   */
  template<typename t_simulation, typename t_interaction_spec>
  class Forcefield : public std::vector<Interaction<t_simulation, 
						    t_interaction_spec> *>
  {
  public:
    /**
     * Constructor
     */
    Forcefield();
    /**
     * Destructor
     */
    ~Forcefield();
    /**
     * calculate all interactions.
     */
    void calculate_interactions(t_simulation &simu);

  protected:

  };
  
} // interaction

// inline functions
#include "forcefield.tcc"

#endif
