/**
 * @file steepest_descent.h
 * steepest descent energy minimisation
 */

#ifndef INCLUDED_STEEPEST_DESCENT_H
#define INCLUDED_STEEPEST_DESCENT_H

namespace algorithm
{
  /**
   * @class Steepest_Descent
   * implements steepest descent energy minimisation
   */
  class Steepest_Descent : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Steepest_Descent() : Algorithm("SteepestDescent") {}

    /**
     * Destructor.
     */
    virtual ~Steepest_Descent(){}
    
    /**
     * steepest descent step.
     */
    virtual int apply(topology::Topology &topo, 
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);

    /**
     * init an algorithm
     * print out input parameter, what it does...
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);
    

  };
  
} // algorithm

#endif

