/**
 * @file algorithm.h
 * base class for algorithms
 */

#ifndef INCLUDED_ALGORITHM_H
#define INCLUDED_ALGORITHM_H

namespace algorithm
{
  /**
   * @class Algorithm
   * base class
   */
  class Algorithm
  {
  public:
    /**
     * Constructor.
     * @param name of the algorithm.
     */
    Algorithm(std::string name) : name(name) {}
    /**
     * Destructor.
     */
    virtual ~Algorithm() {};
    
    /**
     * apply the algorithm
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim) {return 0;};
    /**
     * name of the algorithm
     */
    std::string name;
  };
}

#endif


