/**
 * @file perturbed_flexible_constraint.h
 * the perturbed flexible constraint algorithm.
 */

#ifndef INCLUDED_PERTURBED_FLEXIBLE_CONSTRAINT_H
#define INCLUDED_PERTURBED_FLEXIBLE_CONSTRAINT_H

namespace algorithm
{
  /**
   * @class Perturbed_Flexible_Constraint
   * implements the flexible constraint algorithm 
   * for perturbed distance constraints.
   */
  template<math::virial_enum do_virial>
  class Perturbed_Flexible_Constraint : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Flexible_Constraint(Flexible_Constraint<do_virial> 
				  const & flexible_constraint);
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Flexible_Constraint();
        
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim);
    
    /**
     * initialize startup positions and velocities
     * if required.
     */
    int init(topology::Topology & topo,
	     configuration::Configuration & conf,
	     simulation::Simulation & sim);

  protected:

    Flexible_Constraint<do_virial> const & m_flexible_constraint;
    
  };
  
} //algorithm

// template methods
#include "perturbed_flexible_constraint.tcc"

#endif
