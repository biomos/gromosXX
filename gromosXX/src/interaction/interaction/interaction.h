/**
 * @file interaction.h
 * the interaction interface.
 */

#ifndef INCLUDED_INTERACTION_H
#define INCLUDED_INTERACTION_H

namespace interaction
{
  /**
   * @class interaction
   * @interface interaction
   * declares the interaction interface.
   */
  template<typename t_simulation>
  class interaction
  {
  public:
    /**
     * calculate the interactions.
     */
    virtual void calculate_interactions(t_simulation &simu) = 0;
  };  
  
} // interaction

#endif
