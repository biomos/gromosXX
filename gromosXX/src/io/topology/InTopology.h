/**
 * @file InTopology.h
 * read in a G96 topology file.
 */

#ifndef INCLUDED_INTOPOLOGY_H
#define INCLUDED_INTOPOLOGY_H

namespace algorithm
{
  template<typename t_simulation>
  class Shake;
}

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
    void read_TOPOLOGY(t_topology &topo);
    /**
     * Read in the harmonic bond parameters.
     */
    template<typename t_simulation, typename t_interaction_spec>
    InTopology & operator>>(interaction::harmonic_bond_interaction<t_simulation, 
			    t_interaction_spec> &hbi);
    /**
     * Read in the bond parameter to Shake.
     */
    template<typename t_simulation>
    InTopology & operator>>(algorithm::Shake<t_simulation> &shake);
    /**
     * Read in the quartic bond parameters.
     */    
    template<typename t_simulation, typename t_interaction_spec>
      InTopology & operator>>(interaction::Quartic_bond_interaction<t_simulation,
			      t_interaction_spec> &qbi);
    
    /**
     * Read in the angle parameters.
     */
    template<typename t_simulation, typename t_interaction_spec>
    InTopology & operator>>(interaction::angle_interaction<t_simulation,
			    t_interaction_spec> &ai);
    
    /**
     * Read in the improper dihedral parameters.
     */
    template<typename t_simulation, typename t_interaction_spec>
    InTopology & operator>>
    (interaction::Improper_dihedral_interaction<t_simulation,
     t_interaction_spec> & ii);
    
    /** 
     * Read in the dihedral parameters.
     */
    template<typename t_simulation, typename t_interaction_spec>
      InTopology & operator>>
      (interaction::Dihedral_interaction<t_simulation,
       t_interaction_spec> & di);
    
    /**
     * Read in the nonbonded interaction types (lennard-jones).
     */
    template<typename t_simulation, typename t_interaction_spec>
    InTopology & operator >>(interaction::
			     Nonbonded_Interaction<t_simulation, t_interaction_spec>
			     &nbi);

  };
  

} // io

// template and inline methods
#include "InTopology.tcc"

#endif
