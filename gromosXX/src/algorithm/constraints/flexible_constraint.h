/**
 * @file flexible_constraint.h
 * the flexible shake algorithm
 */

#ifndef INCLUDED_FLEXIBLE_CONSTRAINT_H
#define INCLUDED_FLEXIBLE_CONSTRAINT_H

namespace algorithm
{
  /**
   * @class Flexible_Constraint
   * calculates the flexible constraint distance
   */
  template<typename math::virial_enum do_virial>
  class Flexible_Constraint : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Flexible_Constraint(double const tolerance = 0.000001,
			int const max_iterations = 1000);

    /**
     * Destructor.
     */
    virtual ~Flexible_Constraint();

    /**
     * initialization.
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim);

    /**
     * apply shake.
     */
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

  protected:
    /**
     * shake tolerance
     */
    double m_tolerance;
    /**
     * max iterations.
     */
    const int m_max_iterations;
    /**
     * bond parameter
     */
    std::vector<interaction::bond_type_struct> m_parameter;

  };
  
  
} // algorithm

#include "flexible_constraint.tcc"

#endif