/**
 * @file perturbed_shake.h
 * the perturbed shake algorithm.
 */

#ifndef INCLUDED_PERTURBED_SHAKE_H
#define INCLUDED_PERTURBED_SHAKE_H

namespace algorithm
{
  /**
   * @class Perturbed_Shake
   * implements the shake algorithm for perturbed distance constraints.
   */
  template<math::virial_enum do_virial>
  class Perturbed_Shake : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Shake(Shake<do_virial> const & shake);
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Shake();
        
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

    Shake<do_virial> const & m_shake;
    
  };
  
} //algorithm

// template methods
#include "perturbed_shake.tcc"

#endif
