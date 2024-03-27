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
 * @file virtual_atom.h
 * 
 */

#ifndef INCLUDED_VIRTUAL_ATOM_H
#define INCLUDED_VIRTUAL_ATOM_H

namespace configuration
{
  class Configuration;
}

namespace topology
{
  class Topology;
}

namespace util
{
/**
   * @enum virtual_type
   * holds the virtual atom type
   */
  enum virtual_type {
    /**
     * 0: explicit atom
     */
    va_explicit = 0,
    /**
     * 1: aliphatic CH1
     */
    va_CH1 = 1,
    /**
     * 2: CH1 (aromatic)
     */
    va_aromatic = 2,
    /**
     * 3: CH2 (non-stereospecific)
     */
    va_CH2 = 3,
    /**
     * 8: TIP4P
     */
    va_tip4p = 8,
    /**
     * 4: CH2 (stereospecific)
     */
    va_stereo_CH2 = 4,
    /**
     * 5: CH3
     */
    va_stereo_CH3 = 5,
    /**
     * 6: CH3 (non-stereospecific, Val, Leu)
     */
    va_CH3 = 6,
    /**
     * 7: 3CH3 (non-stereospecific)
     */
    va_3CH3 = 7,
    /**
     * -1: centre of geometry
     */
    va_cog = -1,
    /**
     * -2: centre of mass
     */
    va_com = -2
  };
  
  /**
   * @class Virtual_Atom
   * get position of a virtual or pseudo atom based on 
   * positions of specified other atoms.
   * distribute forces on virtual atom over the atoms used
   * to calculate the virtual site.
   * 
   */

  class Virtual_Atom
  {
    public:
    /**
     * Constructor
     */
    Virtual_Atom();
    /**
     * Constructor
     */
    Virtual_Atom(virtual_type type, std::vector<int> atom, double dish = 0.1, 
		 double disc = 0.153,  int orientation = 0);
    /**
     * calculate the position of the virtual atom
     */
    math::Vec pos(configuration::Configuration & conf, topology::Topology & topo)const;
    /**
     * distribute force f of virtual atom on the real atoms.
     */
    void force(configuration::Configuration & conf, topology::Topology & topo, math::Vec const f)const;
    /**
     * distribute force f of virtual atom on the real atoms (for eds).
     */
    void force(configuration::Configuration & conf, topology::Topology  & topo, math::Vec const f, math::VArray & force)const;
    /**
     * real atom accessor
     */
    int atom(int i)const { assert(i >= 0 && i < int(m_atom.size())); return m_atom[i]; } 
    /**
     * number of atoms that define virtual atom
     */
    int size()const { return m_atom.size(); }
    /**
     * accessor to the type
     */
    virtual_type type() const { return m_type; }
    /**
     * accessor to the charge
     */
    double charge() const { return m_charge; }
    /**
     * accessor to the IAC
     */
    int iac() const { return m_iac; }
    /**
     * set charge
     */
    void set_charge(double charge) { m_charge = charge; }
    /**
     * set iac
     */
    void set_iac(int iac) { m_iac = iac; }

  private:
    /**
     * calculate the position of the virtual site
     */
    template<math::boundary_enum B>
    void _pos(math::VArray const & position, topology::Topology const & topo, math::Box const & box, 
	      math::Vec & p)const;
    /**
     * distribute the force of the virtual site on the real atoms
     */
    template<math::boundary_enum B>
    void _force(math::VArray const & position, topology::Topology const & topo,  math::Box const & box, 
		math::Vec const & f, math::VArray & force)const;
    /**
     * type of the virtual atom
     */
    virtual_type m_type;
    /**
     * atoms to specify the virtual site
     */
    std::vector<int> m_atom;
    /**
     * C-H bond length
     */
    double m_dish;
    /**
     * C-C bond length
     */
    double m_disc;
    /**
     * orientation
     */
    int m_orientation;
    /**
     * charge
     */
    double m_charge;
    /**
     * IAC
     */
    int m_iac;
  };


    /**
     * @class Virtual_Atoms_Group
     * holds the information about Virtual atoms in topology.
     */
  class Virtual_Atoms_Group
  {
  public:
      
    /**
     * virtual atoms accessor.
     */
    std::map<unsigned int, Virtual_Atom> & atoms() {return m_atom;}

    /**
     * const virtual atoms accessor.
     */
    std::map<unsigned int, Virtual_Atom> const & atoms()const {return m_atom;}

    /**
     * virtual atom accessor
     */
    Virtual_Atom & atom(unsigned int i) {return m_atom[i];}
      
  private:
      
    /**
     * the virtual atoms.
     */
    std::map<unsigned int, Virtual_Atom> m_atom; 

  };
}    

#endif
