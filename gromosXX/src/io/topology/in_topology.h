/**
 * @file in_topology.h
 * read in a G96 topology file.
 */

#ifndef INCLUDED_IN_TOPOLOGY_H
#define INCLUDED_IN_TOPOLOGY_H


#include "../ifp.h"
#include "../instream.h"

namespace io {

  /**
   * @class In_Topology
   * reads in a topology file and parses
   * it into topology::Topology
   * @sa topology::Topology
   */
  class In_Topology : public GInStream, public IFP {

  public:
    /**
     * Default constructor.
     */
    In_Topology() {}
    /**
     * destructor
     */
    virtual ~In_Topology() {}
    
    /**
     * Constructor.
     */
    In_Topology(std::istream& is) : GInStream(is) { readStream(); }

    /**
     * set the stream
     */
    void stream(std::istream& is) { GInStream::stream(is); readStream(); }

    /**
     * Read in a G96 topology into the topology.
     */
    void read(topology::Topology &topo, simulation::Parameter &param,
	      std::ostream & os = std::cout);

    /**
     * Read in the harmonic bond parameter.
     */
    virtual void read_harmonic_bonds(std::vector<interaction::bond_type_struct> &b,
				     std::ostream & os = std::cout);

    /**
     * Read in the quartic bond parameter.
     */    
    virtual void read_g96_bonds(std::vector<interaction::bond_type_struct> &b,
				std::ostream & os = std::cout);
    
    /**
     * Read in the bond angle parameter.
     */    
    virtual void read_angles(std::vector<interaction::angle_type_struct> &a,
			     std::ostream & os = std::cout);

    /**
     * Read in the harmonic bond angle parameter.
     */    
    virtual void read_harm_angles(std::vector<interaction::angle_type_struct> &a,
				  std::ostream & os = std::cout);

    /**
     * Read in the improper dihedral parameter.
     */
    virtual void read_improper_dihedrals(std::vector<interaction::improper_dihedral_type_struct> &i,
					 std::ostream & os = std::cout);

    /**
     * Read in the dihedral parameter.
     */
    virtual void read_dihedrals(std::vector<interaction::dihedral_type_struct> &d,
				std::ostream & os = std::cout);
    
    /**
     * Read in the nonbonded interaction types (lennard-jones).
     */
    virtual void read_lj_parameter(std::vector<std::vector
				   <interaction::lj_parameter_struct> > 
				   & lj_parameter,
				   std::ostream & os = std::cout);

    /**
     * Read in the nonbonded interaction types (lennard-jones).
     */
    virtual void read_cg_parameter(std::vector<std::vector
				   <interaction::lj_parameter_struct> > 
				   & cg_parameter,
				   std::ostream & os = std::cout);

    /**
     * Read in the nonbonded interaction types (SASA).
     */
    virtual void read_sasa_parameter(topology::Topology & topo, std::vector
                                    <topology::sasa_parameter_struct>
                                    & sasa_parameter);
    
    /**
     * length of strings allowed
     */
    static const unsigned int MAX_NAME = 5;
    
  private:
    /**
     * solute bond types
     */
    int num_solute_bondtypes;

  };

} // io

#endif
