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
   * @class Leap_Frog
   * implements the leap frog algorithm.
   */
  template<typename t_simulation,
	   typename t_thermostat = algorithm::Berendsen_Thermostat>
  class Leap_Frog
  {
  public:
    typedef t_simulation simulation_type;
    
    /**
     * Constructor.
     */
    Leap_Frog();
    
    /**
     * Leap frog step.
     */
    template<typename t_interaction_spec>
    void step(t_simulation &sim,
	      interaction::Forcefield<t_simulation, t_interaction_spec> &ff, 
	      double const dt);

  protected:
    /**
     * Leap frog velocities.
     */
    void velocities(typename t_simulation::system_type &sys,
		    typename t_simulation::topology_type &topo,
		    double const dt);
    
    /**
     * Leap frog positions.
     */
    void positions(typename t_simulation::system_type &sys,
		   double const dt);
    
    /**
     * The thermostat.
     */
    t_thermostat m_thermostat;

  };
  
} // algorithm

// template definition
#include "leap_frog.tcc"

#endif

