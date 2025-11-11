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

// gcore_Dihedral.h

#ifndef INCLUDED_GCORE_DIHEDRAL
#define INCLUDED_GCORE_DIHEDRAL

namespace gcore{

  /**
   * Class Dihedral
   * Purpose: contains a gromos96 dihedral angle specification
   *
   * Description:
   * Constaint the atoms and type forming a dihedral angle. Atoms are 
   * sorted in such a way that the two middle atoms are in ascending order
   *
   * @class Dihedral
   * @author R. Buergi
   * @ingroup gcore
   * @sa gcore::DihedralType
   * @sa gcore::MoleculeTopology
   */
class Dihedral{
  int d_a[4];
  int d_type;
  // not implemented
  Dihedral();
 public:
  /**
   * Dihedral constructor. The atoms are stored such that c<d.
   * @param a,b,c,d atom numbers defining a dihedral angle
   */
  Dihedral(int a, int b, int c, int d, bool warn = true);
  /**
   * Dihedral copy constructor.
   * @param & Dihedral to be copied
   */
  Dihedral(const Dihedral &);
  /**
   * Dihedral deconstructor
   */
  ~Dihedral(){}
  /**
   * Member operator = copies one Dihedral to the other
   */
  Dihedral &operator=(const Dihedral &);
  /**
   * Method to set the type of the Dihedral
   */
  void setType(int i){d_type = i;}
  /**
   * Accessor, returns the i-th atom in the dihedral as a const (0<i<3)
   */
  int operator[](int i)const{return d_a[i];}
  /**
   * Accessor, returns the i-th atom in the dihedral (0<i<3)
   */
  int &operator[](int i){return d_a[i];}
  /**
   * Accessor, returns the type of the Dihedral
   */
  int type()const{return d_type;}
};
/**
 * @relates Dihedral
 * Operator < compares two dihedral angles a and b
 *
 * a is defined to be smaller than b if <ol>
 * <li> The atom number of the second atom in Dihedral a is smaller than
 *      the atom number of the second atom in Dihedral b (a[1]<b[1]), or
 * <li> The second atoms in a and b have the same atom number, but the 
 *      third atom number in a is lower than the third atom number in b
 *      (a[1]==b[1] && a[2]<b[2]), or
 * <li> The middle atom numbers of a and b are all the same, but the first
 *      atom number in a is lower than the first atom number in b
 *      (a[1]==b[1] && a[2]==b[2] && a[0]<b[0]), or
 * <li> The first three atom numbers in a and b are all the same, but the
 *      fourth atom number in a is lower than the fourth atom number in b
 *      (a[0]==b[0] && a[1]==b[1] && a[2]==b[2] && a[3]<b[3]), or  
 * <li> All atom numbers in a and b are the same, but a has a lower bond
 *      type than b (a[0]==b[0] && a[1]==b[1] && a[2]==b[2] && a[3]==b[3] 
 *      && a.type()<b.type())
 *</ol>
 * @param a,b the Dihedrals to be compared
 * @return 1 if a<b 0 otherwise
 */
int operator<(const Dihedral &a, const Dihedral &b);

}
#endif
