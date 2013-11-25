/**
 * @file perturbed_distance_restraint_interaction.h
 * perturbed distance restraining
 */

#ifndef INCLUDED_PERTURBED_DISTANCE_RESTRAINT_INTERACTION_H
#define INCLUDED_PERTURBED_DISTANCE_RESTRAINT_INTERACTION_H

namespace interaction
{
  /**
   * @class Perturbed_Distance_Restraint_Interaction
   * calculates the perturbed distance restraining interaction
   */
  class Perturbed_Distance_Restraint_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Distance_Restraint_Interaction() : Interaction("PerturbedDistanceRestraint"),
              exponential_term(0.0) {}
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Distance_Restraint_Interaction() {
      DEBUG(2, "Perturbed_Distance_Restraint_Interaction: destructor");
    }

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
    double exponential_term;
    
  };
  
} // interaction

#endif