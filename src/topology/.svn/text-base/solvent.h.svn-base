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
	: name(n), residue_nr(r), iac(i), mass(m), charge(c),
          polarisability(0.0), coscharge(0.0),
          damping_level(0.0), damping_power(0.0),
          gamma(0.0), gamma_j(1), gamma_k(1) {};
      
      std::string name;
      int residue_nr;
      int iac;
      double mass;
      double charge;
      double polarisability;
      double coscharge;
      double damping_level;
      double damping_power;
      double gamma;
      int gamma_j;
      int gamma_k;
    };
    
    /**
     * accessor to the solvent atom information.
     */
    std::vector<atom_struct> & atoms(){return m_atom;}

    /**
     * const accessor to the solvent atom information.
     */
    const std::vector<atom_struct> & atoms() const {return m_atom;}

    /**
     * accessor to a single atom of solvent atom struct.
     */
    atom_struct & atom(unsigned int i){assert(i < m_atom.size()); return m_atom[i];}
    /**
     * const accessor to a single atom of solvent atom struct.
     */
    atom_struct const & atom(unsigned int i)const{assert(i < m_atom.size()); return m_atom[i];}

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

  
    
