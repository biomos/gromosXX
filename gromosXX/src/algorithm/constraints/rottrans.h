/**
 * @file rottrans.h
 * roto-translational constraints
 */

#ifndef INCLUDED_ROTTRANS_H
#define INCLUDED_ROTTRANS_H

namespace algorithm
{
  /**
   * @class Rottrans_Constraints
   * implements roto-translational constraints
   */
  class Rottrans_Constraints : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Rottrans_Constraints(std::string const name = "Rottrans_Constraints");

    /**
     * Destructor.
     */
    virtual ~Rottrans_Constraints();
    
    /**
     * apply roto-translational constraints
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim);

    /**
     * initialise
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);

  protected:

  };
  
} //algorithm

#endif


