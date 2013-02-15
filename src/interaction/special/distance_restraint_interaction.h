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
    Distance_Restraint_Interaction() : Interaction("DistanceRestraint"), 
        exponential_term(0.0) {}
    
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
		     bool quiet = false);
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);

  protected:
    /**
     * exponential term used in averaging
     */
    double exponential_term;
    
  };
  
} // interaction

#endif
