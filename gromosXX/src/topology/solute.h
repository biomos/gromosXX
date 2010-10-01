/**
 * @file solute.h
 * the solute class.
 */

#ifndef INCLUDED_SOLUTE_H
#define INCLUDED_SOLUTE_H

namespace topology
{
  /**
   * @class Solute
   * holds Solute information.
   */
  class Solute : public Compound
  {
  public:
    /**
     * @struct atom_struct
     * holds the information.
     */
    struct atom_struct
    {
      atom_struct(std::string n, int r) : name(n),residue_nr(r){};
      
      std::string name;
      int residue_nr;
    };

    /**
     * Constructor
     */
    Solute() : Compound(){};

    /**
     * accessor to all atom information.
     */
    std::vector<atom_struct> & atoms(){return m_atom;}

    /**
     * accessor to the atom information.
     */
    atom_struct & atom(unsigned int i){assert (i < m_atom.size()); return m_atom[i];}
    /**
     * const accessor to the atom information.
     */
    atom_struct const & atom(unsigned int i)const{assert (i < m_atom.size()); return m_atom[i];}
    
    /**
     * add a soluteatom at the end of the list.
     */
    void add_atom(std::string name, int residue_nr){
      m_atom.push_back(atom_struct(name, residue_nr));
      ++m_num_atoms;
    }

    /**
     * bond accessor.
     */
    std::vector<two_body_term_struct> & bonds(){return m_bond;}

    /**
     * const bond accessor.
     */
    std::vector<two_body_term_struct> const & bonds()const
    {
      return m_bond;
    }

    /**
     * coarse grained bond accessor.
     */
    std::vector<two_body_term_struct> & cgbonds(){return m_cgbond;}

    /**
     * const coarse grained bond accessor.
     */
    std::vector<two_body_term_struct> const & cgbonds()const
    {
      return m_cgbond;
    }
    
    /**
     * angle accessor.
     */
    std::vector<three_body_term_struct> & angles(){return m_angle;}
    /**
     * const angle accessor.
     */
    std::vector<three_body_term_struct> const & angles()const{return m_angle;}
    
    /**
     * improper dihedral accessor.
     */
    std::vector<four_body_term_struct> & improper_dihedrals(){
      return m_improper_dihedral;}
    /**
     * const improper dihedral accessor.
     */
    std::vector<four_body_term_struct> const & improper_dihedrals()const{
      return m_improper_dihedral;}

    /**
     * dihedral accessor.
     */
    std::vector<four_body_term_struct> & dihedrals(){
      return m_dihedral;}
    /**
     * const dihedral accessor.
     */
    std::vector<four_body_term_struct> const & dihedrals()const{
      return m_dihedral;}

    /**
     * crossdihedral accessor.
     */
    std::vector<eight_body_term_struct> & crossdihedrals(){
      return m_crossdihedral;}
    /**
     * const crossdihedral accessor.
     */
    std::vector<eight_body_term_struct> const & crossdihedrals()const{
      return m_crossdihedral;}

  private:

    /**
     * atom information.
     */
    std::vector<atom_struct> m_atom;
    
    /**
     * the bonds.
     */
    std::vector<two_body_term_struct> m_bond;

    /**
     * the coarse grained bonds.
     */
    std::vector<two_body_term_struct> m_cgbond;

    /**
     * the angles.
     */
    std::vector<three_body_term_struct> m_angle;
    
    /**
     * the improper dihedral.
     */
    std::vector<four_body_term_struct> m_improper_dihedral;
    
    /**
     * the dihedral
     */
    std::vector<four_body_term_struct> m_dihedral;

    /**
     * the crossdihedral
     */
    std::vector<eight_body_term_struct> m_crossdihedral;
    
  };
  
} // topology

#endif
