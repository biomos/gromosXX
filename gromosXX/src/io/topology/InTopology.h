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
   * it into simulation::Topology
   * @sa simulation::Topology
   */
  class InTopology : public GInStream {

  public:
    /**
     * Default constructor.
     */
    InTopology() {}
    /**
     * Constructor.
     */
    InTopology(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a G96 topology into the topology.
     */
    template<typename t_topology>
    InTopology & operator>>(t_topology &topo);
    /**
     * Read in the harmonic bond parameters.
     */
    template<typename t_simulation>
    InTopology & operator>>(interaction::harmonic_bond_interaction<t_simulation> &hbi);
    /**
     * Read in the quartic bond parameters.
     */
    template<typename t_simulation>
      InTopology & operator>>(interaction::Quartic_bond_interaction<t_simulation> &qbi);
    
    /**
     * Read in the angle parameters.
     */
    template<typename t_simulation>
    InTopology & operator>>(interaction::angle_interaction<t_simulation> &ai);
    
    /**
     * Read in the improper dihedral parameters.
     */
    template<typename t_simulation>
    InTopology & operator>>
    (interaction::Improper_dihedral_interaction<t_simulation> & ii);
    
    /** 
     * Read in the dihedral parameters.
     */
    template<typename t_simulation>
      InTopology & operator>>
      (interaction::Dihedral_interaction<t_simulation> & di);
    
    /**
     * Read in the nonbonded interaction types (lennard-jones).
     */
    template<typename t_simulation, typename t_pairlist, typename t_innerloop>
    InTopology & operator >>(interaction::
			     Nonbonded_Interaction<t_simulation, t_pairlist, t_innerloop>
			     &nbi);

  };
  

} // io

// template and inline methods
#include "InTopology.tcc"

#endif
