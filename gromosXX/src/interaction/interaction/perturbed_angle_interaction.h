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
  template<typename t_simulation>
  class Perturbed_Angle_Interaction : public Interaction<t_simulation>
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Angle_Interaction(angle_interaction<t_simulation> 
				       &angle_interaction);
    /**
     * Destructor.
     */
    virtual ~Perturbed_Angle_Interaction();
    /**
     * calculate the interactions (force and energy, lambda derivative)
     */
    virtual void calculate_interactions(t_simulation &sim);
    
  protected:
    angle_interaction<t_simulation> &m_angle_interaction;
  };
  
} // interaction

// template methods
#include "perturbed_angle_interaction.tcc"

#endif
