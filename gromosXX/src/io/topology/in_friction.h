/**
 * @file in_posres.h
 * read in a position restraining file.
 */

#ifndef INCLUDED_IN_FRICTION_H
#define INCLUDED_IN_FRICTION_H

#include <gromosXX/io/instream.h>

namespace io {

  /**
   * @class In_Friction
   * reads in a position restraining file
   */
  class In_Friction : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Friction() {}
    /**
     * Constructor.
     */
    In_Friction(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a friction specification file.
     */
    void read(topology::Topology &topo,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);

  };

} // io

#endif
