/**
 * @file perturbed_angle_restraint_interaction.h
 * perturbed angle restraining
 */

#ifndef INCLUDED_PERTURBED_ANGLE_RESTRAINT_INTERACTION_H
#define INCLUDED_PERTURBED_ANGLE_RESTRAINT_INTERACTION_H

namespace interaction
{
  /**
   * @class Perturbed_Angle_Restraint_Interaction
   * calculates the perturbed angle restraining interaction
   */
  class Perturbed_Angle_Restraint_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Angle_Restraint_Interaction() : Interaction("PerturbedAngleRestraint") {}
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Angle_Restraint_Interaction() {}

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
