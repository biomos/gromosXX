/**
 * @file perturbed_angle_interaction.h
 * perturbed angle interaction.
 */

#ifndef INCLUDED_PERTURBED_ANGLE_INTERACTION
#define INCLUDED_PERTURBED_ANGLE_INTERACTION

namespace interaction
{
  /**
   * @class Perturbed_Angle_Interaction
   * calculates the perturbed angle interactions.
   */
  template<typename t_interaction_spec>
  class Perturbed_Angle_Interaction : 
    public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Angle_Interaction(Angle_Interaction<t_interaction_spec> 
				&angle_interaction)
      : Interaction("PerturbedAngle"),
	m_interaction(angle_interaction)
    {}
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Angle_Interaction() {}

    /**
     * calculate the interactions (force and energy, lambda derivative)
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
    
  protected:
    Angle_Interaction<t_interaction_spec> &m_interaction;
  };
  
} // interaction

// template methods
#include "perturbed_angle_interaction.cc"

#endif
