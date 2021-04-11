/**
 * @file propagat_forces.h
 * Propagates forces of the virtual atoms
 */

#ifndef PROPVIRT_H
#define	PROPVIRT_H

namespace algorithm
{
   /**
   * @class Prepare_VirtualAtoms
   * calculates total energies, updates the averages
   */
  class Propagate_Forces : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Propagate_Forces() : Algorithm("PropagateForces"){}

    /**
     * Destructor.
     */
    virtual ~Propagate_Forces(){}
    
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
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false){ 
                 os << "PROPAGATE FORCES OF VIRTUAL ATOMS\nEND\n";
                 return 0; }
  };
  
}
#endif