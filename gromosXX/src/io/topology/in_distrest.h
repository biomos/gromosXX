/**
 * @file in_distrest.h
 * read in a distrance restraints file.
 */

#ifndef INCLUDED_IN_DISTREST_H
#define INCLUDED_IN_DISTREST_H

#include <gromosXX/io/instream.h>

namespace io {

  /**
   * @class In_Distrest
   * reads in a position restraints file
   */
  class In_Distrest : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Distrest() {}
    /**
     * Constructor.
     */
    In_Distrest(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a position restraints file.
     */
    void read(topology::Topology &topo,
	      configuration::Configuration & conf,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);

  };

} // io

#endif
