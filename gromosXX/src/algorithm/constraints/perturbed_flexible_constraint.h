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
  class Perturbed_Flexible_Constraint : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Flexible_Constraint(Flexible_Constraint const & flexible_constraint);
    
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
	     simulation::Simulation & sim,
	     bool quiet = false);

  protected:

    Flexible_Constraint const & m_flexible_constraint;
    
  };
  
} //algorithm

#endif
