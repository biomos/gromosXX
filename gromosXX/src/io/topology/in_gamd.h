/**
 * @file in_gamd.h
 * read in a gamd atoms file.
 */
/**
 * @page GAMD acceleration groups format
 * @date 19-11-2020
 *
 * A Gaussian acceleration MD specification file may contain the following
 * blocks:
 * - @ref title
 * - @ref gamdatoms
 */

#ifndef INCLUDED_IN_GAMD_H
#define INCLUDED_IN_GAMD_H

#include "../instream.h"

namespace io {

  /**
   * @class In_GAMD
   * reads in a GAMD atoms specification file
   */
  class In_GAMD : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_GAMD() {}
    /**
     * Constructor.
     */
    In_GAMD(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a xray restraining specification file.
     */
    void read(topology::Topology &topo,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);

  };
  
} // io

#endif