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
  class Position_Restraint_Interaction : public Interaction
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
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false) 
    {
      os << "Position restraint interaction\n";
      return 0;
    };
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);

  protected:
    
  };
  
} // interaction

#endif
