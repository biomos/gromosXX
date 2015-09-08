/**
 * @file in_rdc.h
 * read in a rdc restraining specification file.
 */
/**
 * @page rdc RDC restraints specification format
 * @date 09-05-2011
 *
 * A RDC restraints specifcation file has to contain the following
 * blocks:
 * - @ref title
 * - @ref conversion
 * - @ref magfieldc
 * - @ref alignt
 * - @ref rdcresspec
 */

#ifndef INCLUDED_IN_RDC_H
#define INCLUDED_IN_RDC_H

#include "../instream.h"

namespace io {

  /**
   * @class In_RDC
   * reads in a RDC restraint specification file.
   */
  class In_RDC : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_RDC() {}
    /**
     * Constructor.
     */
    In_RDC(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a RDC restraining file.
     */
    void read(topology::Topology &topo,
	      configuration::Configuration & conf,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);

  };

} // io

#endif

