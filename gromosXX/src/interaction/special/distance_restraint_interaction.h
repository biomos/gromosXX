/**
 * @file distance_restraint_interaction.h
 * distance restraining
 */

#ifndef INCLUDED_DISTANCE_RESTRAINT_INTERACTION_H
#define INCLUDED_DISTANCE_RESTRAINT_INTERACTION_H

namespace interaction
{
  /**
   * @class distance_restraint_interaction
   * calculates the distance restraining interaction
   */
  class Distance_Restraint_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Distance_Restraint_Interaction() : Interaction("DistanceRestraint") {}
    
    /**
     * Destructor.
     */
    virtual ~Distance_Restraint_Interaction() {}

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false) 
    {
      if (!quiet)
	os << "Distance restraint interaction\n";
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
