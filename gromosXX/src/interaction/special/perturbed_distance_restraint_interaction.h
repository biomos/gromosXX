/**
 * @file perturbed_distance_restraint_interaction.h
 * perturbed distance restraining
 */

#ifndef INCLUDED_PERTURBED_DISTANCE_RESTRAINT_INTERACTION_H
#define INCLUDED_PERTURBED_DISTANCE_RESTRAINT_INTERACTION_H

namespace interaction
{
  /**
   * @class perturbed_distance_restraint_interaction
   * calculates the distance restraining interaction
   */
  class Perturbed_Distance_Restraint_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Distance_Restraint_Interaction() : Interaction("PerturbedDistanceRestraint") {}
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Distance_Restraint_Interaction() {}

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false) 
    {
      os << "Perturbed distance restraint interaction\n";
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
