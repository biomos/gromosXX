/**
 * @file perturbed_solute.h
 * the perturbed part of the solute.
 */

#ifndef INCLUDED_PERTURBED_SOLUTE_H
#define INCLUDED_PERTURBED_SOLUTE_H

namespace topology
{
  /**
   * @class Perturbed_Solute
   * holds the information about perturbation of
   * the solute.
   */
  class Perturbed_Solute
  {
  public:
    /**
     * const perturbed bonds.
     */
    std::vector<perturbed_two_body_term_struct> const & bonds()const {return m_bond;}

    /**
     * perturbed bonds.
     */
    std::vector<perturbed_two_body_term_struct> & bonds() {return m_bond;}

    /**
     * const perturbed coarse grained bonds.
     */
    std::vector<perturbed_two_body_term_struct> const & cgbonds()const {return m_cgbond;}

    /**
     * perturbed coarse grained bonds.
     */
    std::vector<perturbed_two_body_term_struct> & cgbonds() {return m_cgbond;}

    /**
     * const perturbed angles.
     */
    std::vector<perturbed_three_body_term_struct> const & angles()const {return m_angle;}
    
    /**
     * perturbed angles.
     */
    std::vector<perturbed_three_body_term_struct> & angles() {return m_angle;}
    
    /**
     * const perturbed improper dihedrals.
     */
    std::vector<perturbed_four_body_term_struct> const & improper_dihedrals()const{
      return m_improper_dihedral;
    }
    
    /**
     * perturbed improper dihedrals.
     */
    std::vector<perturbed_four_body_term_struct> & improper_dihedrals(){
      return m_improper_dihedral;
    }
    
    /**
     * const perturbed dihedrals.
     */
    std::vector<perturbed_four_body_term_struct> const & dihedrals()const
    {
      return m_dihedral;
    }
    
    /**
     * perturbed dihedrals.
     */
    std::vector<perturbed_four_body_term_struct> & dihedrals()
    {
      return m_dihedral;
    }
    
    /**
     * perturbed atoms accessor.
     */
    std::map<unsigned int, Perturbed_Atom> & atoms() {return m_atom;}

    /**
     * const perturbed atoms accessor.
     */
    std::map<unsigned int, Perturbed_Atom> const & atoms()const {return m_atom;}

    /**
     * perturbed atom accessor
     */
    Perturbed_Atom & atom(unsigned int i) {return m_atom[i];}

    /**
     * perturbed atompairs.
     */
    std::vector<perturbed_two_body_term_struct> & atompairs(){return m_atompair;}
    
    /**
     * const perturbed atompairs.
     */
    std::vector<perturbed_two_body_term_struct> const & atompairs()const {return m_atompair;}

    /**
     * perturbed distance constraints accessor.
     */
    std::vector<perturbed_two_body_term_struct> & distance_constraints() {
      return m_distance_constraint;
    }

    /**
     * perturbed distance constraints const accessor.
     */
    std::vector<perturbed_two_body_term_struct> const & distance_constraints()const{
      return m_distance_constraint;
    }
    
  private:
    /**
     * the perturbed bonds.
     */
    std::vector<perturbed_two_body_term_struct> m_bond;

    /**
     * the perturbed coarse grained bonds.
     */
    std::vector<perturbed_two_body_term_struct> m_cgbond;

    /**
     * the perturbed angles.
     */
    std::vector<perturbed_three_body_term_struct> m_angle;

    /**
     * the perturbed improper dihedrals.
     */
    std::vector<perturbed_four_body_term_struct> m_improper_dihedral;
    
    /**
     * the perturbed dihedrals.
     */
    std::vector<perturbed_four_body_term_struct> m_dihedral;
    
    /**
     * the perturbed atoms.
     */
    std::map<unsigned int, Perturbed_Atom> m_atom;
    
    /**
     * the perturbed atompairs.
     */
    std::vector<perturbed_two_body_term_struct> m_atompair;
    
    /**
     * the perturbed distance constraints.
     */
    std::vector<perturbed_two_body_term_struct> m_distance_constraint;

  };
      

  /**
   * @class EDS_Perturbed_Solute
   * holds the information about EDS-perturbation of
   * the solute.
   */
  class EDS_Perturbed_Solute
  {
  public:
     
    /**
     * perturbed atoms accessor.
     */
    std::map<unsigned int, EDS_Perturbed_Atom> & atoms() {return m_atom;}

    /**
     * const perturbed atoms accessor.
     */
    std::map<unsigned int, EDS_Perturbed_Atom> const & atoms()const {return m_atom;}

    /**
     * perturbed atom accessor
     */
    EDS_Perturbed_Atom & atom(unsigned int i) {return m_atom[i];}
      
  private:
     
    /**
     * the perturbed atoms.
     */
    std::map<unsigned int, EDS_Perturbed_Atom> m_atom; 
 
  };
      

} // topology


#endif
