/**
 * @file forcefield.h
 * the Forcefield class.
 */

#ifndef INCLUDED_FORCEFIELD_H
#define INCLUDED_FORCEFIELD_H

namespace configuration{
	class Configuration;
}
namespace topology{
	class Topology;
}
namespace simulation{
	class Simulation;
}

/**
 * @namespace interaction
 * namespace that contains the classes to
 * handle the interactions between the particles.
 * (energies, forces).
 */
namespace interaction
{
  /**
   * @class Forcefield
   * contains the specific interactions.
   * (clear does not call them (i guess) -- sorry, don't know what this means anymore)
   * Strategy Pattern.
   */
  class Forcefield : public std::vector<Interaction *>,
		     public algorithm::Algorithm
  {
  public:
    /**
     * Constructor
     */
    Forcefield() 
      : std::vector<Interaction *>(),
	algorithm::Algorithm("Forcefield") {}
    /**
     * Destructor
     */
    ~Forcefield();
    /**
     * initialise
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);

    /**
     * calculate all interactions.
     */
    int calculate_interactions(topology::Topology & topo,
			       configuration::Configuration & conf,
			       simulation::Simulation & sim);


    /**
     * let the forcefield be used as an algorithm
     */
    int apply(topology::Topology & topo,
	      configuration::Configuration & conf,
	      simulation::Simulation & sim)
    {
      return calculate_interactions(topo, conf, sim);
    }
    
    virtual void print_timing(std::ostream & os);

    /**
     * const interaction accessor.
     */
    Interaction const * interaction(std::string name)const;

    /**
     * interaction accessor.
     */
    Interaction * interaction(std::string name);
    
  protected:

  };
  
} // interaction

#endif
