/**
 * @file eds_distance_restraint_interaction.h
 * eds perturbed distance restraining
 */

#ifndef INCLUDED_EDS_DISTANCE_RESTRAINT_INTERACTION_H
#define INCLUDED_EDS_DISTANCE_RESTRAINT_INTERACTION_H

namespace interaction
{
  /**
   * @class eds_distance_restraint_interaction
   * calculates the eds perturbed distance restraining interaction
   */
  class Eds_Distance_Restraint_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Eds_Distance_Restraint_Interaction() : Interaction("EdsDistanceRestraint"), 
        exponential_term(0.0) {}
    
    /**
     * Destructor.
     */
    virtual ~Eds_Distance_Restraint_Interaction() {}

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
