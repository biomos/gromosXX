/**
 * @file in_distanceres.h
 * read in a distrance restraints file.
 */
/**
 * @page disres distance restraints format
 * @date 28-10-2008
 *
 * A distance restraints specifcation file may contain the following blocks:
 * - @ref title
 * - @ref distanceresspec
 * - @ref pertdisresspec
 * - @ref mdisresspec
 * - @ref dfresspec
 * - @ref pertdfresspec
 */

#ifndef INCLUDED_IN_DISTANCERES_H
#define INCLUDED_IN_DISTANCERES_H

#include "../instream.h"

namespace io {

  /**
   * @class In_Distanceres
   * reads in a position restraints file
   */
  class In_Distanceres : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Distanceres() {}
    /**
     * Constructor.
     */
    In_Distanceres(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a distance restraints file.
     */
    void read(topology::Topology &topo,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);
    /**
     * Maximum number of atoms that can be specified to define a virtual atom
     */
    static const unsigned int MAX_ATOMS = 4;
  };

} // io

#endif
