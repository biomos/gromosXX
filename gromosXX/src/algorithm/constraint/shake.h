/**
 * @file shake.h
 * the shake algorithm.
 */

#ifndef INCLUDED_SHAKE_H
#define INCLUDED_SHAKE_H

namespace algorithm
{
  /**
   * @class shake
   * implements the shake algorithm.
   */
  template<typename t_simulation>
  class shake
  {
  public:
    typedef t_simulation simulation_type;
    /**
     * Constructor.
     */
    shake();

    /**
     * shake solvent.
     */
    int solvent(typename simulation_type::topology_type &topo,
		typename simulation_type::system_type &sys,
		double dt);

  private:
    double tolerance;
    const int max_iterations;
  };
  
} //algorithm

// template methods
#include "shake.tcc"

#endif
  
