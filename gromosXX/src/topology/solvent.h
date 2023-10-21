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

  
    
