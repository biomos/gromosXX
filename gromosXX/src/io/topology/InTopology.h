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
    InTopology(std::istream& is) : GInStream(is) { read_stream(); };
    /**
     * Read in a G96 topology into the topology.
     */
    InTopology & operator>>(simulation::topology &topo);
    /**
     * Read in the harmonic bond parameters.
     */
    template<typename t_simulation>
    InTopology & operator>>(interaction::harmonic_bond_interaction<t_simulation> &hbi);
    /**
     * Read in the nonbonded interaction types (lennard-jones).
     */
    template<typename t_simulation>
    InTopology & operator >>(interaction::nonbonded_interaction<t_simulation> &nbi);

  private:
    /**
     * read the entire stream and store the blocks in the map.
     */
    void read_stream();

    std::map<std::string, std::vector<std::string> > m_block;
    
  };
  

} // io

// template and inline methods
#include "InTopology.tcc"

#endif
