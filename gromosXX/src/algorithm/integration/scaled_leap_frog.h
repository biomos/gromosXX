/**
 * @file scaled_leap_frog.h
 * the leap frog algorithm.
 */

#ifndef INCLUDED_SCALED_LEAP_FROG_H
#define INCLUDED_SCALED_LEAP_FROG_H

namespace algorithm
{
  /**
   * @class Scaled_Leap_Frog_Velocity
   * implements the leap frog algorithm for the velocities.
   */
  class Scaled_Leap_Frog_Velocity : public Leap_Frog_Velocity
  {
  public:
/**
     * Constructor.
     */
    Scaled_Leap_Frog_Velocity() :
    Leap_Frog_Velocity() { name = "Scaled_Leap_Frog_Velocity"; }
    /**
     * Destructor.
     */
    virtual ~Scaled_Leap_Frog_Velocity(){};
    
    /**
     * Scaled leap frog step.
     */
    virtual int apply(topology::Topology &topo, 
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);


  };

  
} // algorithm

#endif

