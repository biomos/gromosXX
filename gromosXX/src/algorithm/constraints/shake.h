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
  template<math::virial_enum do_virial>
  class Shake : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Shake(double const tolerance = 0.000001, int const max_iterations = 1000);
    /**
     * Destructor.
     */
    virtual ~Shake();
        
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim);
    
    /**
     * set the tolerance.
     */
    void tolerance(double const tol);

    /**
     * tolerance.
     */
    double const tolerance()const {return m_tolerance;}
    /**
     * max iterations.
     */
    int const max_iterations()const {return m_max_iterations;}

    /**
     * the const bond type parameter.
     */
    std::vector<interaction::bond_type_struct> const &parameter()const
    {
      return m_parameter;
    }
    /**
     * the bond type parameter.
     */
    std::vector<interaction::bond_type_struct> & parameter()
    {
      return m_parameter;
    }

    /**
     * initialize startup positions and velocities
     * if required.
     */
    int init(topology::Topology & topo,
	     configuration::Configuration & conf,
	     simulation::Simulation & sim);

  protected:

    double m_tolerance;
    const int m_max_iterations;
    
    std::vector<interaction::bond_type_struct> m_parameter;

  };
  
} //algorithm

// template methods
#include "shake.tcc"

#endif
