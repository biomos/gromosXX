/**
 * @file forcefield.h
 * the Forcefield class.
 */

#ifndef INCLUDED_FORCEFIELD_H
#define INCLUDED_FORCEFIELD_H

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
   * @TODO are the destructors called? properly?
   * clear does not call them (i guess).
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
    

  protected:

  };
  
} // interaction

#endif
