/**
 * @file position_restraint_interaction.h
 * position restraining
 */

#ifndef INCLUDED_POSITION_RESTRAINT_INTERACTION_H
#define INCLUDED_POSITION_RESTRAINT_INTERACTION_H

namespace interaction
{
  /**
   * @class position_restraint_interaction
   * calculates the position restraining interaction
   */
  template<typename t_interaction_spec>
  class Position_Restraint_Interaction : 
    public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Position_Restraint_Interaction() : Interaction("PositionRestraint") {}
    
    /**
     * Destructor.
     */
    virtual ~Position_Restraint_Interaction() {}

    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);

  protected:
    
  };
  
} // interaction

// template methods
#include "position_restraint_interaction.tcc"

#endif
