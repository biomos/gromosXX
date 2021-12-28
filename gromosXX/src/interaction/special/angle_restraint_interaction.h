/**
 * @file angle_restraint_interaction.h
 * angle restraining
 */

#ifndef INCLUDED_ANGLE_RESTRAINT_INTERACTION_H
#define INCLUDED_ANGLE_RESTRAINT_INTERACTION_H

namespace interaction
{
  /**
   * @class angle_restraint_interaction
   * calculates the angle restraining interaction
   */
  class Angle_Restraint_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Angle_Restraint_Interaction() : Interaction("AngleRestraint") {}
    
    /**
     * Destructor.
     */
    virtual ~Angle_Restraint_Interaction() {}

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false);

    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
  };
  
} // interaction

#endif
