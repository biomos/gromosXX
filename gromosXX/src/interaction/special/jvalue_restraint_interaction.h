/**
 * @file jvalue_restraint_interaction.h
 * J-value restraining
 */

#ifndef INCLUDED_JVALUE_RESTRAINT_INTERACTION_H
#define INCLUDED_JVALUE_RESTRAINT_INTERACTION_H

namespace interaction
{
  /**
   * @class Jvalue_Restraint_Interaction
   * calculates the J-Value restraining interaction
   */
  template<typename t_interaction_spec>
  class Jvalue_Restraint_Interaction : 
    public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Jvalue_Restraint_Interaction() : Interaction("JValueRestraint") {}
    
    /**
     * Destructor.
     */
    virtual ~Jvalue_Restraint_Interaction() {}

    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
    
  };
  
} // interaction

// template methods
#include "jvalue_restraint_interaction.cc"

#endif
