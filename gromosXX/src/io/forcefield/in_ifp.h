// gio_InParameter.h

#ifndef INCLUDED_IN_IFP_H
#define INCLUDED_IN_IFP_H

#include "../ifp.h"
#include "../instream.h"

namespace io{

  /**
   * Class In_IFP
   * defines an instream that can read a GROMOS96 ifp-file
   *
   * The GROMOS96 ifp file is read in and stored in a GromosForceField
   * format. This means that vdw-parameters are already calculated to 
   * the individual pairs, taking all special tables in the manual into 
   * account'
   *
   * @class In_IFP
   * @author B.C. Oostenbrink
   */
  class In_IFP : public GInStream, public IFP {
  public:
    /**
     * Constructor
     */
    In_IFP(std::istream& is) : GInStream(is) { readStream(); }
    /**
     * Destructor
     */
    virtual ~In_IFP() {}

    //////////////////////////////////////////////////
    // INTERFACE
    //////////////////////////////////////////////////

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
    
  private:
    
  };
} // io

#endif
