/**
 * @file in_friction.h
 * read in a friction specification file.
 */

/**
 * @page friction friction specification format
 * @date 28-10-2008
 *
 * A friction specifcation file may contain the following
 * blocks:
 * - @ref title
 * - @ref frictionspec
 */

#ifndef INCLUDED_IN_FRICTION_H
#define INCLUDED_IN_FRICTION_H

#include "../instream.h"

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
