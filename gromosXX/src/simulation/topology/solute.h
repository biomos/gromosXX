/**
 * @file solute.h
 * the solute class.
 */

#ifndef INCLUDED_SOLUTE_H
#define INCLUDED_SOLUTE_H

namespace simulation
{
  /**
   * @class solute
   * holds solute information.
   */
  class Solute : public compound
  {
  public:
    /**
     * @struct soluteatom_struct
     * holds the information.
     */
    struct atom_struct
    {
      std::string name;
      int residue_nr;
    };

    /**
     * Constructor
     */
    Solute();

    /**
     * accessor to all atom information.
     */
    std::vector<atom_struct> & atoms();

    /**
     * accessor to the atom information.
     */
    atom_struct & atom(size_t i);
    
    /**
     * add a soluteatom at the end of the list.
     */
    void add_atom(std::string name, int residue_nr);

    /**
     * bond accessor.
     */
    std::vector<Bond> & bonds();

    /**
     * add all bonds to the solute constraint vector and
     * remove them from the bond vector.
     */
    void add_bond_length_constraints(std::vector<interaction::bond_type_struct> const &param);
    
    /**
     * add bonds connecting an atom of type iac to the
     * constraint vector and remove from the bond vector.
     */
    void add_bond_length_constraints(int iac, std::vector<int> const &atom_iac,
				    std::vector<interaction::bond_type_struct> const &param);
    
    /**
     * add bonds connecting an atom of mass mass to the
     * constraint vector and remove from the bond vector.
     */
    void add_bond_length_constraints(double mass, math::SArray const &atom_mass,
				    std::vector<interaction::bond_type_struct> const & param);

    /**
     * angle accessor.
     */
    std::vector<Angle> & angles();
    
    /**
     * improper dihedral accessor.
     */
    std::vector<Improper_Dihedral> & improper_dihedrals();
    /**
     * dihedral accessor.
     */
    std::vector<Dihedral> & dihedrals();
    
    
  private:
    /**
     * atom information.
     */
    std::vector<atom_struct> m_atom;
    
    /**
     * the bonds.
     */
    std::vector<Bond> m_bond;
    
    /**
     * the angles.
     */
    std::vector<Angle> m_angle;
    
    /**
     * the improper dihedral.
     */
    std::vector<Improper_Dihedral> m_improper_dihedral;
    
    /**
     * the dihedral
     */
    std::vector<Dihedral> m_dihedral;
    
  };
  
} // simulation

// inline methods
#include "solute.tcc"

#endif
