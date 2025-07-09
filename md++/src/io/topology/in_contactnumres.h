/**
 * @file in_contactnumres.h
 * read in a collective variable restraints file.
 */
/**
 * @page colvarres colvar restraints format
 * @date 28-10-2008
 *
 * A colvar restraints specification file may contain the following blocks:
 * - @ref title
 * - @ref in_contactnumresspec
 */

#ifndef INCLUDED_IN_CONTACTNUMRES_H
#define INCLUDED_IN_CONTACTNUMRES_H

#include "../instream.h"

namespace io {

  /**
   * @class In_Colvarres
   * reads in a collective variable restraints file
   */
  class In_Contactnumres : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Contactnumres() {}
    /**
     * Constructor.
     */
    In_Contactnumres(std::istream& is) : GInStream(is) { readStream(); };
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
