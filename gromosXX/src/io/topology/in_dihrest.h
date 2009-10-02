/**
 * @file in_dihrest.h
 * read in a dihedral restraints file.
 */
/**
 * @page dihrest dihedral restraints format
 * @date 28-10-2008
 *
 * A dihedral restraints/constraints specifcation file may contain the following
 * blocks:
 * - @ref title
 * - @ref dihedralresspec
 * - @ref pertdihresspec
 */
#ifndef INCLUDED_IN_DIHREST_H
#define INCLUDED_IN_DIHREST_H

#include "../instream.h"

namespace io {

  /**
   * @class In_Dihrest
   * reads in a dihedral restraints file
   */
  class In_Dihrest : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Dihrest() {}
    /**
     * Constructor.
     */
    In_Dihrest(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a position restraints file.
     */
    void read(topology::Topology &topo,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);
  };
} // io

#endif
