/**
 * @file src/interaction/interaction/interaction.h
 * the interaction interface.
 */

#ifndef INCLUDED_INTERACTION_H
#define INCLUDED_INTERACTION_H

namespace interaction
{
  /**
   * @class Interaction
   * @interface Interaction
   * declares the interaction interface.
   */
  template<typename t_simulation>
  class Interaction
  {
  public:
    /**
     * Constructor.
     */
    Interaction(std::string name) : name(name) {};
    /**
     * Destructor.
     */
    virtual ~Interaction(){};
    /**
     * the name of the interaction.
     * can be used to identify a special class.
     */
    std::string name;
    /**
     * calculate the interactions.
     */
    virtual void calculate_interactions(t_simulation &simu) = 0;
  };  
  
} // interaction

#endif
