/**
 * @file algorithm.h
 * the Algorithm class.
 */

#ifndef INCLUDED_ALGORITHM_ALGORITHM_H
#define INCLUDED_ALGORITHM_ALGORITHM_H

/**
 * @namespace algorithm
 * namespace that contains the classes to
 * handle the interactions between the particles.
 * (energies, forces).
 */
namespace algorithm
{
  /**
   * @class Algorithm_Sequence
   * contains the specific algorithms.
   */
  class Algorithm_Sequence : public std::vector<Algorithm *>
  {
  public:
    /**
     * Constructor
     */
    Algorithm_Sequence();
    /**
     * Destructor
     */
    ~Algorithm_Sequence();
    /**
     * calculate all interactions.
     */
    int run(topology::Topology &topo, 
	    configuration::Configuration &conf,
	    simulation::Simulation &sim);

  protected:

  };
  
} // algorithm


#endif
