/**
 * @file in_dihrest.h
 * read in a dihedral restraints file.
 */

#ifndef INCLUDED_IN_DIHREST_H
#define INCLUDED_IN_DIHREST_H

#include <gromosXX/io/instream.h>

namespace io {

  /**
   * @class In_Dihrest
   * reads in a dihedral restraints file
   */
  class In_Dihrest : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Dihrest() {}
    /**
     * Constructor.
     */
    In_Dihrest(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a position restraints file.
     */
    void read(topology::Topology &topo,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);
  };
} // io

#endif
