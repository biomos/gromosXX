/**
 * @file runge_kutta.h
 * the runge kutta algorithm.
 */

#ifndef INCLUDED_RUNGE_KUTTA_H
#define INCLUDED_RUNGE_KUTTA_H

namespace algorithm
{
  /**
   * @class runge_kutta
   * implements the runge-kutta integration scheme.
   */
  template<typename t_simulation>
  class runge_kutta
  {
  public:
    typedef t_simulation simulation_type;
    /**
     * 4th order Runge-Kutta integration step.
     */
    void step(t_simulation &sim,
	      interaction::forcefield<t_simulation> &ff,
	      double const dt);
    
  protected:
    /**
     * Store the intermediate function evaluations (and gradients)
     * for the positions.
     */
    math::VArray dx1, dx2, dx3, dx4;
    /**
     * Store the intermediate function evaluations (and gradients)
     * for the velocity.
     */
    math::VArray dv1, dv2, dv3, dv4;
  };
  
} // algorithm

// template definitions
#include "runge_kutta.tcc"

#endif

