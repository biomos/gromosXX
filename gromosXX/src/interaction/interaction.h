/**
 * @file src/interaction/interaction.h
 * the interaction interface.
 */

#ifndef INCLUDED_INTERACTION_H
#define INCLUDED_INTERACTION_H

namespace configuration{
	class Configuration;
}
namespace topology{
	class Topology;
}
namespace simulation{
	class Simulation;
}

namespace interaction
{
  /**
   * @class Interaction
   * @interface Interaction
   * declares the interaction interface.
   */
  class Interaction
  {
  public:
    /**
     * Constructor.
     */
    Interaction(std::string name) : name(name), m_timing(0.0) {};
    /**
     * Destructor.
     */
    virtual ~Interaction(){};
    /**
     * the name of the interaction.
     * can be used to identify a special class.
     */
    std::string name;
    /**
     * initialise
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false) = 0;
    // { return 0; }
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim) = 0;

    /**
     * timing information.
     */
    virtual void print_timing(std::ostream & os)
    {
      os << "        "
	 << std::setw(36) << std::left << name
	 << std::setw(20) << m_timing << "\n";
    }

  protected:
    /**
     * store time used in algorithm.
     */
    double m_timing;

  };  
  
} // interaction

#endif
