/**
 * @file leap_frog.h
 * the leap frog algorithm.
 */

#ifndef INCLUDED_LEAP_FROG_H
#define INCLUDED_LEAP_FROG_H

/**
 * @namespace algorithm
 * namespace that contains the implementations of
 * the various algorithms.
 */
namespace algorithm
{
  /**
   * @class leap_frog
   * implements the leap frog algorithm.
   */
  template<typename t_simulation>
  class leap_frog
  {
  public:
    typedef t_simulation simulation_type;
    
    /**
     * Leap frog step.
     */
    static void step(t_simulation &sim,
		     interaction::Forcefield<t_simulation> &ff, 
		     double const dt);

  protected:
    /**
     * Leap frog velocities.
     */
    static void velocities(typename t_simulation::system_type &sys,
			   typename t_simulation::topology_type &topo,
			   double const dt);
    
    /**
     * Leap frog positions.
     */
    static void positions(typename t_simulation::system_type &sys,
			  double const dt);

  };
  
} // algorithm

// template definition
#include "leap_frog.tcc"

#endif

