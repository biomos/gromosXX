/**
 * @file angle_interaction.h
 * angle interaction.
 */

#ifndef INCLUDED_ANGLE_INTERACTION_H
#define INCLUDED_ANGLE_INTERACTION_H

namespace interaction
{
  /**
   * @class angle_interaction
   * calculates the angle interactions.
   */
  template<typename t_simulation>
  class angle_interaction : public Interaction<t_simulation>
  {
  public:
    /**
     * Constructor.
     */
    angle_interaction();
    /**
     * Destructor.
     */
    virtual ~angle_interaction();
    /**
     * calculate the interactions.
     */
    virtual void calculate_interactions(t_simulation &simu);
    /**
     * add bondangle type.
     */
    void add(angle_type_struct s);
    /**
     * add bond type.
     */
    void add(double K, double cos0);
    /**
     * the angle type parameters.
     */
    std::vector<angle_type_struct> const & parameter()const;
    
  protected:
    std::vector<angle_type_struct> m_angle_parameter;
    
  };
  
} // interaction

// template methods
#include "angle_interaction.tcc"

#endif
