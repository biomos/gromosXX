/**
 * @file in_zaxisoribias.h
 * read in a z-axis orientation bias file.
 */
/**
 * @page zaxisoribias z-axis orientation bias format
 * @date 18-06-2019
 *
 * A z-axis orientation bias specification file may contain the following blocks:
 * - @ref title
 * - @ref zaxisoribiasspec
 */

#ifndef INCLUDED_IN_ZAXISORIBIAS_H
#define INCLUDED_IN_ZAXISORIBIAS_H

#include "../instream.h"

namespace io {

  /**
   * @class In_Zaxisoribias
   * reads in a z-axis orientation biass file
   */
  class In_Zaxisoribias : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Zaxisoribias() {}
    /**
     * Constructor.
     */
    In_Zaxisoribias(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a z-axis orientation restraints file.
     */
    void read(topology::Topology &topo,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);

    /**
     * read z-axis orientation bias specification block.
     */
    void read_ZAXISORIBIASSPEC(topology::Topology &topo, simulation::Simulation &sim, std::ostream & os = std::cout);

    /**
     * Maximum number of atoms that can be specified to define a virtual atom
     */
    static const unsigned int MAX_ATOMS = 4;
  };

} // io

#endif
