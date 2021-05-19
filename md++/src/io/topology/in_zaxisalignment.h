/**
 * @file in_zaxisalignment.h
 * read in a z-axis angle restraints file.
 */
/**
 * @page zangleres z-axis angle restraints format
 * @date 18-06-2019
 *
 * A z-axis angle restraints specification file may contain the following blocks:
 * - @ref title
 * - @ref zalignmentresspec
 */

#ifndef INCLUDED_IN_ZALIGNMENTRES_H
#define INCLUDED_IN_ZALIGNMENTRES_H

#include "../instream.h"

namespace io {

  /**
   * @class In_Zaxisalignment
   * reads in a z-axis angle restraints file
   */
  class In_Zaxisalignment : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Zaxisalignment() {}
    /**
     * Constructor.
     */
    In_Zaxisalignment(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a z-axis angle restraints file.
     */
    void read(topology::Topology &topo,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);

    /**
     * read z-axis angle restraint specification block.
     */
    void read_ZALIGNMENTRESSPEC(topology::Topology &topo, simulation::Simulation &sim, std::ostream & os = std::cout);

    /**
     * Maximum number of atoms that can be specified to define a virtual atom
     */
    static const unsigned int MAX_ATOMS = 4;
  };

} // io

#endif
