/**
 * @file in_leus.h
 * read in a local elevation umbrella samping file.
 */
/**
 * @page leus LE-US format
 * @date 04-022-2009
 *
 * A local elevation file may contain the following
 * blocks:
 * - @ref title
 * - @ref localelevspec
 *
 * @sa leusbias
 *
 * @page leusdb LE-US database format
 * @date 04-022-2009
 *
 * A local elevation umbrella samping database file may contain the following
 * blocks:
 * - @ref title
 * - @ref leusbias
 *
 * @sa leus
 */

#ifndef INCLUDED_IN_LEUS_H
#define INCLUDED_IN_LEUS_H

#include "../instream.h"

namespace io {

  /**
   * @class In_Localelevspec
   * reads in a local elevation specification file
   */
  class In_Localelevspec : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Localelevspec() {}
    /**
     * Constructor.
     */
    In_Localelevspec(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a local elevation specification file.
     */
    void read(topology::Topology &topo,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);
  };

  /**
   * @class In_LEUSBias
   * reads in a local elevation umbrella samping database file
   */
  class In_LEUSBias : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_LEUSBias() {}
    /**
     * Constructor.
     */
    In_LEUSBias(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a local elevation umbrella sampling database file.
     */
    void read(topology::Topology &topo,
              configuration::Configuration & conf,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);
  };

} // io

#endif


