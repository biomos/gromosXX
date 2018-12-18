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
	 * Read topology blocks
	*/
	void read_block_TYPE(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_PHYSICALCONSTANTS(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_RESNAME(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_ATOMTYPENAME(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_SOLUTEATOM(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_SOLUTEPOLARISATION(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_CGSOLUTE(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_LJEXCEPTIONS(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_BONDH(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_BOND(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_BONDDP(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_CONSTRAINT(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_BONDANGLEH(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_BONDANGLE(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_IMPDIHEDRAL(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_IMPDIHEDRALH(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_DIHEDRAL(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_DIHEDRALH(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_CROSSDIHEDRAL(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_CROSSDIHEDRALH(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_VIRTUALGRAIN(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_SOLUTEMOLECULES(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_TEMPERATUREGROUPS(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_PRESSUREGROUPS(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);

    void read_SOLVENT_blocks(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_SOLVENTPOLARISATION(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os, topology::Solvent &s);
    void read_block_SOLVENTCONSTR(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os, topology::Solvent &s);

    void read_block_BONDSTRETCHTYPE(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_BONDTYPE(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_HARMBONDTYPE(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);

    void read_block_BONDANGLEBENDTYPE(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_BONDANGLETYPE(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_HARMBONDANGLETYPE(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);

    void read_block_DIHEDRALTYPE(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_TORSDIHEDRALTYPE(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);
    void read_block_IMPDIHEDRALTYPE(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os);

    /**
     * Read in the bond parameters.
     */
    void read_bond_types(topology::Topology& topo,
       simulation::Parameter &param,
       std::ostream & os);

    /**
    * Read in the bond angle parameters.
    */
    void read_bondangle_types(topology::Topology& topo,
      simulation::Parameter &param,
      std::ostream & os);

   /**
   * Read in the dihedral angle parameters.
   */
   void read_dihedral_types(topology::Topology& topo,
     simulation::Parameter &param,
     std::ostream & os);

      /**
      * Read in the charged virtual sites.
      */
      virtual void read_offsite_chg( std::vector
                                     <interaction::off_site_struct >
                                     & offsite_parameter,topology::Topology & topo);
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
    /**
     * solute angle types
     */
    int num_solute_angletypes;
    /**
     * solute angle types
     */
    int num_solute_dihedraltypes;
    /**
     * solute angle types
     */
    int num_solute_impropertypes;

  };

} // io

#endif
