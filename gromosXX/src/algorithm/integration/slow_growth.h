/**
 * @file slow_growth.h
 * slow growth algorithm
 */

#ifndef INCLUDED_SLOW_GROWTH_H
#define INCLUDED_SLOW_GROWTH_H

namespace algorithm
{
  /**
   * @class Slow_Growth
   * implements slow growth.
   * mainly updates the topology for a changed
   * lambda value.
   */
  class Slow_Growth : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Slow_Growth() : Algorithm("SlowGrowth") {}

    /**
     * Destructor.
     */
    virtual ~Slow_Growth() {}
    
    /**
     * do the lambda change.
     */
    virtual int apply(topology::Topology &topo, 
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);

  };
  
} // algorithm

#endif
