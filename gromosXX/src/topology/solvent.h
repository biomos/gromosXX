/**
 * @file solvent.h
 * solvent information.
 * allows multiple solvents!
 */

#ifndef INCLUDED_SOLVENT_H
#define INCLUDED_SOLVENT_H

namespace topology
{
  /**
   * @class Solvent
   * holds solvent information
   * (for one solute).
   */
  class Solvent : public Compound
  {
  public:
    /**
     * @struct atom_struct
     * the solvent atom information.
     */
    struct atom_struct
    {
      atom_struct(std::string n, int r, int i, double m, double c)
	: name(n), residue_nr(r), iac(i), mass(m), charge(c) {};
      
      std::string name;
      int residue_nr;
      int iac;
      double mass;
      double charge;
    };
    
    /**
     * accessor to the solvent atom information.
     */
    std::vector<atom_struct> & atoms(){return m_atom;}

    /**
     * accessor to a single atom of solvent atom struct.
     */
    atom_struct & atom(size_t i){assert(i < m_atom.size()); return m_atom[i];}
    /**
     * const accessor to a single atom of solvent atom struct.
     */
    atom_struct const & atom(size_t i)const{assert(i < m_atom.size()); return m_atom[i];}

    /**
     * add a solventatom.
     */
    void add_atom(std::string name, int res_nr, int iac, 
		  double mass, double charge){
      m_atom.push_back(atom_struct(name, res_nr, iac, mass, charge));
      ++m_num_atoms;
    }

  private:
    std::vector<atom_struct> m_atom;

  };
  
} // topology

#endif

  
    
