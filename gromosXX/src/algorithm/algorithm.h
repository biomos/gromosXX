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
    Algorithm(std::string name) : name(name), m_timing(0.0) {}

    /**
     * Destructor.
     */
    virtual ~Algorithm() {};
    
    /**
     * init an algorithm
     * print out input parameter, what it does...
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim) { return 0; }
    
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

    /**
     * timing information.
     */
    virtual void print_timing(std::ostream & os)
    {
      os << std::setw(40) << std::left << name
	 << std::setw(20) << m_timing << "\n";
    }
    
  protected:
    /**
     * store time used in algorithm.
     */
    double m_timing;
    
  };
}

#endif


