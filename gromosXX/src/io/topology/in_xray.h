/**
 * @file in_xray.h
 * read in a xray restraining file.
 */
/**
 * @page xrayres xray restraints format
 * @date 03-02-2009
 *
 * A xray restraints specifcation file may contain the following
 * blocks:
 * - @ref title
 * - @ref xrayresspec
 * - @ref xrayrespara
 * - @ref elementmap
 */

#ifndef INCLUDED_IN_XRAYRES_H
#define INCLUDED_IN_XRAYRES_H

#include "../instream.h"

namespace io {

  /**
   * @class In_Xrayresspec
   * reads in a xray restraining specification file
   */
  class In_Xrayresspec : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Xrayresspec() {}
    /**
     * Constructor.
     */
    In_Xrayresspec(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a xray restraining specification file.
     */
    void read(topology::Topology &topo,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);

  };
  
} // io

#endif


