/**
 * @file in_angrest.h
 * read in a angle restraints file.
 */
/**
 * @page angrest angle restraints format
 * @date 24-05-2019
 *
 * A angle restraints/constraints specification file may contain the following
 * blocks:
 * - @ref title
 * - @ref angresspec
 * - @ref pertangresspec
 */
#ifndef INCLUDED_IN_ANGREST_H
#define INCLUDED_IN_ANGREST_H

#include "../instream.h"

namespace io {

  /**
   * @class In_Angrest
   * reads in a angle restraints file
   */
  class In_Angrest : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Angrest() {}
    /**
     * Constructor.
     */
    In_Angrest(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a position restraints file.
     */
    void read(topology::Topology &topo,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);

    void read_ANGRESSPEC(topology::Topology &topo,
        									simulation::Simulation &sim,
        									std::ostream & os);
    void read_PERTANGRESSPEC(topology::Topology &topo,
									simulation::Simulation &sim,
									std::ostream & os);
  };
} // io

#endif
