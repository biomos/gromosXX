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
    bond & bonds();

    /**
     * add all bonds to the solute constraint vector and
     * remove them from the bond vector.
     */
    void add_bond_length_constraints();
    
    /**
     * add bonds connecting an atom of type iac to the
     * constraint vector and remove from the bond vector.
     */
    void add_bond_length_constraint(int iac);
    
    /**
     * add bonds connecting an atom of mass mass to the
     * constraint vector and remove from the bond vector.
     */
    void add_bond_length_constraint(double mass);

  private:
    /**
     * atom information.
     */
    std::vector<atom_struct> m_atom;
    
    /**
     * the bonds.
     */
    bond m_bond;
    
  };
  
} // simulation

// inline methods
#include "solute.tcc"

#endif
