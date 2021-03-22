/**
 * @file analyze.h
 * re-evaluate a trajectory
 */

#ifndef INCLUDED_ANALYZE_H
#define INCLUDED_ANALYZE_H

namespace algorithm
{
  /**
   * @class Analyze_Step
   * implements analyzation of a trajectory
   */
  class Analyze_Step : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Analyze_Step(std::string trajectory);

    /**
     * Destructor.
     */
    virtual ~Analyze_Step(){}
    
    /**
     * analyze step.
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
    
  private:
    /**
     * trajectory to analyze
     */
    io::In_Configuration m_trajectory;

    /**
     * trajectory file
     */
   io::igzstream m_trajectory_file;
    
  };
  
   /**
   * @class Analyze_Exchange
   * implements analyzation of a trajectory. Does nothing but exchanging the
   * states. Substitutes leap-frog in re-evaluation.
   */
  class Analyze_Exchange : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Analyze_Exchange();

    /**
     * Destructor.
     */
    ~Analyze_Exchange(){}
    
    /**
     * Do exchange states.
     */
    virtual int apply(topology::Topology &topo, 
		      configuration::Configuration &conf,
                      simulation::Simulation &sim);

    /**
     * Initialization
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
                     bool quiet = false);
  };
  
} // algorithm

#endif
