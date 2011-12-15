/**
 * @file   in_bsleus.h
 * Read in the file defining the B&S-LEUS Umbrella
 */

#ifndef IN_BSLEUS_H
#define	IN_BSLEUS_H

#include "../instream.h"

namespace io {
  /**
   * @class In_BSLEUS reads in the definition of the B&S-LEUS Umbrella 
   * definition.
   */
  class In_BSLEUS : public GInStream {
  public:
    /**
     * Constructor
     * @param is the file to read in.
     */
    In_BSLEUS (std::istream &is) : GInStream(is) {readStream();}
    /**
     * Read in the definition file of the B&S-LEUS algorithm
     * @param topo the topology
     * @param sim the simulation parameter
     * @param os  
     */
    void read(topology::Topology &topo,
          configuration::Configuration &conf,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);
  };
}
#endif	/* IN_BSLEUS_H */

