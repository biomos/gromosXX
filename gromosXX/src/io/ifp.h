/**
 * @file ifp.h
 * interaction function parameter read-in interface
 */

#ifndef INCLUDED_IFP_H
#define INCLUDED_IFP_H

namespace io
{
  /**
   * @class IFP
   * interface to read in interaction function parameters.
   */
  class IFP
  {
  public:

    /**
     * Read in the harmonic bond parameter.
     */
    virtual void read_harmonic_bonds(std::vector<interaction::bond_type_struct> &b) = 0;
    /**
     * Read in the quartic bond parameter.
     */    
     virtual void read_g96_bonds(std::vector<interaction::bond_type_struct> &b) = 0;
    /**
     * Read in the bond angle parameter.
     */    
    virtual void read_angles(std::vector<interaction::angle_type_struct> &a) = 0;

    /**
     * Read in the improper dihedral parameter.
     */
    virtual void read_improper_dihedrals(std::vector<interaction::improper_dihedral_type_struct> &i) = 0;

    /**
     * Read in the dihedral parameter.
     */
    virtual void read_dihedrals(std::vector<interaction::dihedral_type_struct> &d) = 0;
    
    /**
     * Read in the nonbonded interaction types (lennard-jones).
     */
    virtual void read_lj_parameter(std::vector<std::vector
				   <interaction::lj_parameter_struct> > 
				   & lj_parameter) = 0;
  };
  
} // namespace io

#endif
