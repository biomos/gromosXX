/**
 * @file position_constraints.h
 * position constraints
 */

#ifndef INCLUDED_POSITION_CONSTRAINTS_H
#define INCLUDED_POSITION_CONSTRAINTS_H

namespace algorithm
{
  /**
   * @class Position_Constraints
   * implements position constraints
   */
  class Position_Constraints : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Position_Constraints(std::string const name = "Position_Constraints");

    /**
     * Destructor.
     */
    virtual ~Position_Constraints();
    
    /**
     * apply position constraints
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


