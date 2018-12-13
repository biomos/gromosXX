/**
 * @file algorithm_sequence.h
 * the Algorithm class.
 */

#ifndef INCLUDED_ALGORITHM_ALGORITHM_H
#define INCLUDED_ALGORITHM_ALGORITHM_H

namespace algorithm
{
  /**
   * @class Algorithm_Sequence
   * contains the specific algorithms.
   * An almost clean implementation of the Strategy pattern.
   */
  class Algorithm_Sequence : public std::vector<Algorithm *>
  {
  public:
    /**
     * Constructor
     */
    Algorithm_Sequence(bool clean = true);

    /**
     * Destructor
     */
    ~Algorithm_Sequence();

    /**
     * init
     */
    int init(topology::Topology &topo, 
	     configuration::Configuration &conf,
	     simulation::Simulation &sim,
	     std::ostream & os = std::cout,
	     bool quiet = false);

    /**
     * calculate all interactions.
     */
    int run(topology::Topology &topo, 
	    configuration::Configuration &conf,
	    simulation::Simulation &sim);

    /**
     * print timing information
     */
    int print_timing(std::ostream & os);

    /**
     * algorithm accessor
     */
    Algorithm * algorithm(std::string name);
    
  protected:
    bool clean;

  };
  
} // algorithm


#endif
