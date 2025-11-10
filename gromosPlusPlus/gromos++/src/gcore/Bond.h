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

// gcore_Bond.h

#ifndef INCLUDED_GCORE_BOND
#define INCLUDED_GCORE_BOND

namespace gcore{

  /**
   * Class Bond
   * Purpose: contains a gromos96 bond
   *
   * Description:
   * Contains the atoms and type making up a bond. The atoms are sorted
   * to have the lowest atom number first.
   *
   * @class Bond
   * @author R. Buergi
   * @ingroup gcore
   * @sa gcore::BondType
   * @sa gcore::MoleculeTopology
   */
class Bond{
  int d_a[2];
  int d_type;
  // not implemented
  Bond();
 public:
  /**
   * Bond constructor
   * constructs a new bond defined by atoms a and b. These atoms are stored
   * such that a<b.
   * @param a,b atom numbers of the atoms making the bond
   */
  Bond(int a, int b, bool warn = true);
  /**
   * Bond copy constructor
   * Constructs a new bond and copies the specied bond into it
   * @param & bond to be copied
   */
  Bond(const Bond &);
  /**
   * Bond deconstructor
   */
  ~Bond(){}
  /**
   * Member operator =, copies one bond into the other
   */
  Bond &operator=(const Bond &b);
  /**
   * Method setType sets the bond type of the bond
   * @param i bond type to be set
   */
  void setType(int i){d_type=i;}
  /**
   * Member operator [], returns the atom number of the i-th atom in the bond
   * @param i atom index in the bond (0 or 1)
   * @return Atom number of the i-th atom in the bond
   */
  int operator[](int i)const{return d_a[i];}
  /**
   * Member operator [], returns the atom number of the i-th atom in the bond
   * @param i atom index in the bond (0 or 1)
   * @return Atom number of the i-th atom in the bond
   */
  int &operator[](int i){return d_a[i];}
  /**
   * Accessor, returns the bond type of the bond
   */
  int type()const{return d_type;}
};
/**
 * @relates Bond
 * Operator < compares two bonds a and b
 *
 * a is defined to be smaller than b if<br><ol>
 * <li> The atom number of the first atom in bond a is lower than the
 *      atom number of the first atom in bond b (a[0]<b[0]), or
 * <li> The atom numbers of the first atoms in both bonds are the same 
 *      but the second atom has a lower atom number in bond a (a[0]==b[0]
 *      && a[1]<b[1]), or
 * <li> All atoms are the same, but the bond type of bond a has a lower 
 *      number than the bond type of bond b (a[0]==b[0] && a[1]==b[1] && 
 *      a.type() < b.type()).
 * </ol>
 * &param a, b bonds to be compared
 * &return 1 if a<b; 0 otherwise
 */
int operator<(const Bond &a, const Bond &b);
}
#endif
