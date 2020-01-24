/**
 * @file in_symrest.h
 * read in a symmetry restraints specification file
 */
/**
 * @page symrest symmetry restraints specification format
 * @date 17-08-2011
 *
 * A symmetry restraints specification file may contain the following
 * blocks:
 * - @ref title
 * - @ref transform
 * - @ref symresspec
 */

#ifndef INCLUDED_IN_SYMREST_H
#define INCLUDED_IN_SYMREST_H

#include "../instream.h"

namespace io {

  /**
   * @class In_Symrest
   * reads in a symmetry restraints specification file
   */
  class In_Symrest : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Symrest() {}
    /**
     * Constructor.
     */
    In_Symrest(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a symmetry restraints specification file.
     */
    void read(topology::Topology &topo,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);

  };
} // io

#endif
