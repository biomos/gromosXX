/**
 * @file in_topology.h
 * read in a G96 topology file.
 */

#ifndef INCLUDED_IN_TOPOLOGY_H
#define INCLUDED_IN_TOPOLOGY_H

namespace io {

  /**
   * @class In_Topology
   * reads in a topology file and parses
   * it into topology::Topology
   * @sa topology::Topology
   */
  class In_Topology : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Topology() {}
    /**
     * Constructor.
     */
    In_Topology(std::istream& is) : GInStream(is) { readStream(); };

    /**
     * set the stream
     */
    void stream(std::istream& is) { GInStream::stream(is); readStream(); }

    /**
     * Read in a G96 topology into the topology.
     */
    void read(topology::Topology &topo, simulation::Parameter &param);

    /**
     * Read in the harmonic bond parameter.
     */
    void read_harmonic_bonds(std::vector<interaction::bond_type_struct> &b);

    /**
     * Read in the quartic bond parameter.
     */    
    void read_g96_bonds(std::vector<interaction::bond_type_struct> &b);

    /**
     * Read in the bond angle parameter.
     */    
    void read_angles(std::vector<interaction::angle_type_struct> &a);

    /**
     * Read in the improper dihedral parameter.
     */
    void read_improper_dihedrals(std::vector<interaction::improper_dihedral_type_struct> &i);

    /**
     * Read in the dihedral parameter.
     */
    void read_dihedrals(std::vector<interaction::dihedral_type_struct> &d);
    
    /**
     * Read in the nonbonded interaction types (lennard-jones).
     */
    void read_lj_parameter(std::vector<std::vector
			   <interaction::lj_parameter_struct> > 
			   & lj_parameter);
  };

} // io

#endif
