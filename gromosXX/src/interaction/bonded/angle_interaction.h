/**
 * @file angle_interaction.h
 * angle interaction.
 */

#ifndef INCLUDED_ANGLE_INTERACTION_H
#define INCLUDED_ANGLE_INTERACTION_H

namespace interaction
{
  /**
   * @class Angle_Interaction
   * calculates the angle interactions.
   */
  template<typename t_interaction_spec>
  class Angle_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Angle_Interaction() : Interaction("Angle") {}
    /**
     * Destructor.
     */
    virtual ~Angle_Interaction() {}
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);

    /**
     * angle interaction parameter.
     */
    std::vector<angle_type_struct> const & parameter()const { return m_parameter; }
    /**
     * angle interaction parameter.
     */
    std::vector<angle_type_struct> & parameter() { return m_parameter; }
    
  protected:
    std::vector<angle_type_struct> m_parameter;
    
  };
  
} // interaction

// template methods
#include "angle_interaction.tcc"

#endif
