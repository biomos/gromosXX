/**
 * @file in_colvarres.h
 * read in a collective variable restraints file.
 */
/**
 * @page colvarres colvar restraints format
 * @date 28-10-2008
 *
 * A colvar restraints specification file may contain the following blocks:
 * - @ref title
 * - @ref contactnumresspec
 */

#ifndef INCLUDED_IN_COLVARRES_H
#define INCLUDED_IN_COLVARRES_H

#include "../instream.h"

namespace io {

  /**
   * @class In_Colvarres
   * reads in a collective variable restraints file
   */
  class In_Colvarres : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Colvarres() {}
    /**
     * Constructor.
     */
    In_Colvarres(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a colvar restraints file.
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
