/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

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
     * const perturbed soft bonds.
     */
    std::vector<perturbed_two_body_term_struct> const & softbonds()const {return m_softbond;}

    /**
     * perturbed soft bonds.
     */
    std::vector<perturbed_two_body_term_struct> & softbonds() {return m_softbond;}
        
    /**
     * const perturbed soft angles.
     */
    std::vector<perturbed_three_body_term_struct> const & softangles()const {return m_softangle;}

    /**
     * perturbed soft angles.
     */
    std::vector<perturbed_three_body_term_struct> & softangles() {return m_softangle;}
        
    /**
     * const perturbed soft improper dihedrals.
     */
    std::vector<perturbed_four_body_term_struct> const & softimpropers()const {return m_softimproper;}

    /**
     * perturbed soft improper dihedrals.
     */
    std::vector<perturbed_four_body_term_struct> & softimpropers() {return m_softimproper;}
    
    /**
     * const perturbed soft bond softness parameter.
     */
    std::vector<double> const & alpha_bond()const {return m_alpha_bond;}

    /**
     * perturbed soft bond softness parameter.
     */
    std::vector<double> & alpha_bond() {return m_alpha_bond;}
    
    /**
     * const perturbed soft angle softness parameter.
     */
    std::vector<double> const & alpha_angle()const {return m_alpha_angle;}

    /**
     * perturbed soft angle softness parameter.
     */
    std::vector<double> & alpha_angle() {return m_alpha_angle;}
    
    /**
     * const perturbed soft improper dihedral softness parameter.
     */
    std::vector<double> const & alpha_improper()const {return m_alpha_improper;}

    /**
     * perturbed soft improper dihedral softness parameter.
     */
    std::vector<double> & alpha_improper() {return m_alpha_improper;}
    
    /**
     * bond_types to be added to the ones in the read topology.
     */
    std::vector<int> const & soft_bond_types()const {return m_soft_bond_types;}
    
    /**
     * bond_types to be added to the ones in the read topology.
     */
    std::vector<int> & soft_bond_types() {return m_soft_bond_types;}
    
    /**
     * angle_types to be added to the ones in the read topology.
     */
    std::vector<int> const & soft_angle_types()const {return m_soft_angle_types;}
    
    /**
     * angle_types to be added to the ones in the read topology.
     */
    std::vector<int> & soft_angle_types() {return m_soft_angle_types;}
    
    /**
     * improper_types to be added to the ones in the read topology.
     */
    std::vector<int> const & soft_improper_types()const {return m_soft_improper_types;}
    
    /**
     * improper_types to be added to the ones in the read topology.
     */
    std::vector<int> & soft_improper_types() {return m_soft_improper_types;}

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
     * the perturbed soft bonds.
     */
    std::vector<perturbed_two_body_term_struct> m_softbond;
    
    /**
     * the perturbed soft angles.
     */
    std::vector<perturbed_three_body_term_struct> m_softangle;
    
    /**
     * the perturbed soft improper dihedrals.
     */
    std::vector<perturbed_four_body_term_struct> m_softimproper;
    
    /**
     * the perturbed soft bond softness parameters.
     */
    std::vector<double> m_alpha_bond;
    
    /**
     * the perturbed soft angle softness parameters.
     */
    std::vector<double> m_alpha_angle;

    /**
     * the perturbed soft improper dihedral softness parameters.
     */
    std::vector<double> m_alpha_improper;
    
    /**
     * bond types to be added for bonds that are made or broken using 
     * the soft potential
     */
    std::vector<int> m_soft_bond_types;
    
    /**
     * angle types to be added for angles that are destroyed
     */
    std::vector<int> m_soft_angle_types;
    
    /**
     * improper types to be added for improper dihedrals that are destroyed
     */
    std::vector<int> m_soft_improper_types;

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
