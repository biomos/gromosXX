/**
 * @file solvent.h
 * solvent information.
 * allows multiple solvents!
 */

#ifndef INCLUDED_SOLVENT_H
#define INCLUDED_SOLVENT_H

namespace simulation
{
  /**
   * @class Solvent
   * holds solvent information
   * (for one solute).
   */
  class Solvent : public compound
  {
  public:
    /**
     * @struct atom_struct
     * the solvent atom information.
     */
    struct atom_struct
    {
      std::string name;
      int residue_nr;
      int iac;
      double mass;
      double charge;
    };
    
    /**
     * accessor to the solvent atom information.
     */
    std::vector<atom_struct> & atoms();

    /**
     * accessor to a single atom of solvent atom struct.
     */
    atom_struct & atom(size_t i);

    /**
     * add a solventatom.
     */
    void add_atom(std::string name, int res_nr, int iac, 
		  double mass, double charge);

  private:
    std::vector<atom_struct> m_atom;

  };
  
} // simulation

// inline methods
#include "solvent.tcc"

#endif

  
    
