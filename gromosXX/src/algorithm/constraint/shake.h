/**
 * @file shake.h
 * the shake algorithm.
 */

#ifndef INCLUDED_SHAKE_H
#define INCLUDED_SHAKE_H

namespace algorithm
{
  /**
   * @class Shake
   * implements the shake algorithm.
   */
  template<typename t_simulation>
  class Shake
  {
  public:
    typedef t_simulation simulation_type;
    /**
     * Constructor.
     */
    Shake(double const tolerance = 0.000001, int const max_iterations = 1000);

    /**
     * shake solute.
     */
    int solute(typename simulation_type::topology_type &topo,
	       typename simulation_type::system_type &sys,
	       double dt);

    /**
     * shake solvent.
     */
    int solvent(typename simulation_type::topology_type &topo,
		typename simulation_type::system_type &sys,
		double dt);

    /**
     * set the tolerance.
     */
    void tolerance(double const tol);

  private:

    bool _shake(typename simulation_type::topology_type const &topo,
		typename simulation_type::system_type &sys,
		int const first, 
		std::vector<bool> &skip_now,
		std::vector<bool> &skip_next,
		std::vector<simulation::compound::distance_constraint_struct>
		& constr);

    double m_tolerance;
    const int max_iterations;
  };
  
} //algorithm

// template methods
#include "shake.tcc"

#endif
  
