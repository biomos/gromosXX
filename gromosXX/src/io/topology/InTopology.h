/**
 * @file InTopology.h
 * read in a G96 topology file.
 */

#ifndef INCLUDED_INTOPOLOGY_H
#define INCLUDED_INTOPOLOGY_H

namespace io {

  /**
   * @class InTopology
   * reads in a topology file and parses
   * it into simulation::topology.
   */
  class InTopology : public GInStream {

  public:
    /**
     * Constructor.
     */
    InTopology(std::istream& is = std::cin) : GInStream(is) {};
    /**
     * Read in a G96 topology into the topology.
     */
    InTopology & operator>>(simulation::topology &topo);
    /**
     * Read in the harmonic bond parameters.
     */
    template<typename t_simulation>
    InTopology & operator>>(interaction::harmonic_bond_interaction<t_simulation> &hbi);

  };
  

} // io

// template and inline methods
#include "InTopology.tcc"

#endif
