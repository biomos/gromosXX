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

// gcore_Constraint.h

#ifndef INCLUDED_GCORE_CONSTRAINT
#define INCLUDED_GCORE_CONSTRAINT

namespace gcore{
  /**
   * Class Constraint
   * Purpose: contains a set of atoms and a distance between them
   *
   * Description:
   * The Constraint contains a set of atoms and a distance between them.
   * It was introduced for the SolventTopology definition
   *
   * @class Constraint
   * @author C. Oostenbrink
   * @ingroup gcore
   * @sa gcore::SolventTopology
   */
class Constraint{
  int d_a[2];
  int d_bondtype;
  double d_dist;
  // not implemented
  Constraint();
 public:
  /**
   * Constraint constructor
   * @param a,b atom numbers of constraint atoms
   */
  Constraint(int a, int b);
  /**
   * Constraint copy constructor
   * @param & Constraint to be copied
   */
  Constraint(const Constraint &);
  /**
   * Constraint deconstructor
   */
  ~Constraint(){}
  /**
   * Member operator = copies one constraint into the other
   */
  Constraint &operator=(const Constraint &b);
  /**
   * Function to set the distance between the atoms
   * @param a the distance to be set
   */
  void setDist(double a){d_dist=a;}

  /**
   * Function to set the bond type corresponding to this constraint
   * (relevant for solute constraints only)
   * @param t the type to be set
   */
  void setType(int t){d_bondtype = t;}

  /**
   * Accessor, returns the atom number of the i-th atom (i=0,1)
   * as a const
   */
  int operator[](int i)const{return d_a[i];}
  /**
   * Accessor, returns the atom number of the i-th atom (i=0,1)
   */
  int &operator[](int i){return d_a[i];}
  /**
   * Accessor, returns the constraint distance
   */
  double dist()const{return d_dist;}

  /**
   * Accessor, returns the bond type corresponding to the constraint
   * (relevant for solute constraints only)
   */
  int bondtype()const{return d_bondtype;}
};

/**
 * @relates Bond
 * Operator < compares two constraints a and b<br>
 * a is defined to be smaller than b if <ol>
 * <li> The atom number of the first atom in constraint a is smaller than
 *      the atom number of the first atom in constraint b (a[0]<b[0]), or
 * <li> The atom numbers of the first atoms in both constraints are the same,
 *      but the second atom has a lower number in constraint a
 *      (a[0]==b[0] && a[1]<b[1])
 * </ol>
 * @param a,b constraints to be compared
 * @return 1 if a<b; 0 otherwise
 */
int operator<(const Constraint &a, const Constraint &b);
}
#endif
