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
 * @file interaction_types.h
 * structures that hold the forcefield parameter.
 */

#ifndef INCLUDED_INTERACTION_TYPES_H
#define INCLUDED_INTERACTION_TYPES_H

namespace interaction
{
  /**
   * @struct interaction_type_struct
   * common base class
   */
  struct interaction_type_struct
  {
  };
  
  /**
   * @struct virtual_atom_type_struct
   * virtual atom types. This stores the distances to compute the virtual atoms.
   */
  struct virtual_atom_type_struct : public interaction_type_struct
  {
    virtual_atom_type_struct(int type, double dis1, double dis2) : type(type), dis1(dis1), dis2(dis2) {}
    int type;
    double dis1;
    double dis2;
  };

  /**
   * @struct bond_type_struct
   * bond types.
   */
  struct bond_type_struct : public interaction_type_struct
  {
    bond_type_struct(double K, double r0) : K(K), r0(r0) {}
    double K;
    double r0;
  };

  /**
   * @struct angle_type_struct
   * angle types.
   */
  struct angle_type_struct : public interaction_type_struct
  {
    angle_type_struct(double K, double cos0) : K(K), cos0(cos0) {}
    double K;
    double cos0;
  };
  
  /**
   * @struct improper_dihedral_type_struct
   * improper dihedral types.
   */
  struct improper_dihedral_type_struct : public interaction_type_struct
  {
    improper_dihedral_type_struct(double K, double q0) : K(K), q0(q0) {}
    double K;
    double q0;
  };
  
  /**
   * @struct dihedral_type_struct
   * dihedral types.
   */
  struct dihedral_type_struct : public interaction_type_struct
  {
    dihedral_type_struct(double K, double cospd, double pd, int m) : K(K), pd(pd), cospd(cospd), m(m) {}
    double K;
    double pd;
    double cospd;
    int m;
  };
  
  /**
   * @struct lj_parameter_struct
   * Lennard Jones interaction parameter.
   */
  struct lj_parameter_struct : public interaction_type_struct
  {
    /**
     * Constructor
     */
    lj_parameter_struct(double c6, double c12, double cs6, double cs12)
      : c6(c6), c12(c12), cs6(cs6), cs12(cs12) {}
    /**
     * Coarse grained constructor (no cs6, cs12)
     */
    lj_parameter_struct(double c6, double c12)
      : c6(c6), c12(c12), cs6(c6), cs12(c12) {}
    /**
     * default constructor
     */
    lj_parameter_struct()
      : c6(0), c12(0), cs6(0), cs12(0) {}
    double c6;
    double c12;
    double cs6;
    double cs12;
  };

}

#endif
