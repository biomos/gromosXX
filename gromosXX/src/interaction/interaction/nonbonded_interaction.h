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
  class Nonbonded_Interaction : 
    public Interaction<t_simulation>,
    public Nonbonded_Base,
    public Nonbonded_Inner_Loop<t_simulation,
				typename t_simulation::system_type>
  {
  public:    
    /**
     * Constructor.
     * @param sim where to store forces and energies
     * (and virial contribution).
     */
    Nonbonded_Interaction(t_simulation &sim);
    
    /**
     * Destructor.
     */
    virtual ~Nonbonded_Interaction();

    /**
     * calculate the interactions.
     */
    virtual void calculate_interactions(t_simulation &sim);

    /**
     * pairlist accessor
     */
    t_pairlist & pairlist();

  protected:
    /**
     * helper class to build the pairlist.
     */
    t_pairlist m_pairlist;
    
    /**
     * calculate the interactions.
     */
    virtual void do_interactions(t_simulation &sim,
				 typename t_pairlist::iterator it, 
				 typename t_pairlist::iterator to);

    /**
     * calculate the 1,4-interactions.
     */
    virtual void do_14_interactions(t_simulation &sim);

    /**
     * calculate the RF contributions for excluded atoms.
     */
    virtual void do_RF_excluded_interactions(t_simulation &sim);
    
  };
  
} // interaction

// template methods
#include "nonbonded_interaction.tcc"

#endif
