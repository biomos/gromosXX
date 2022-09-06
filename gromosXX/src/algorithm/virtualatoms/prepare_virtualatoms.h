/**
 * @file prepare_virtualatoms.h
 * Prepares the virtual atoms for the simulation
 */

#ifndef PREPVIRT_H
#define	PREPVIRT_H

namespace interaction
{
  class Forcefield;
}
namespace algorithm
{
   /**
   * @class Prepare_VirtualAtoms
   * calculates total energies, updates the averages
   */
  class Prepare_VirtualAtoms : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Prepare_VirtualAtoms(): Algorithm("PrepareVirtualAtoms") {}

    /**
     * Destructor.
     */
    virtual ~Prepare_VirtualAtoms(){}
    
    /**
     * calculate new positions of the virtual atoms with nonbonded interactions
     */
    virtual int apply(topology::Topology &topo, 
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);

    /**
     * init
     * extend force and position arrays to hold information of 
     * virtual atoms with nonbonded parameters
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		      std::ostream &os = std::cout, bool quiet = false){
           return 0;
         }
  };
  
}
#endif