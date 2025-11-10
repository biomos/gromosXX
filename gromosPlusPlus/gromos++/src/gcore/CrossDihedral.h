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

// gcore_CrossDihedral.h

#ifndef INCLUDED_GCORE_CROSSDIHEDRAL
#define INCLUDED_GCORE_CROSSDIHEDRAL

namespace gcore{

  /**
   * Class CrossDihedral
   * Purpose: contains a cross dihedral angle specification
   *
   * Description:
   * Contains the atoms and type forming a cross dihedral angle. Atoms are
   * sorted in such a way that the two middle atoms are in ascending order and
   * that the first dihedral is smaller than the second.
   *
   * @class CrossDihedral
   * @author N.Schmid
   * @ingroup gcore
   * @sa gcore::Dihedral
   * @sa gcore::DihedralType
   * @sa gcore::MoleculeTopology
   */
class CrossDihedral{
  int d_a[8];
  int d_type;
  // not implemented
  CrossDihedral();
 public:
  /**
   * CrossDihedral constructor. The atoms are stored such that c<d and g<h.
   * @param a,b,c,d,e,f,g,h atom numbers defining a cross dihedral
   */
  CrossDihedral(int a, int b, int c, int d, int e, int f, int g, int h);
  /**
   * CrossDihedral copy constructor.
   * @param & CrossDihedral to be copied
   */
  CrossDihedral(const CrossDihedral &);
  /**
   * CrossDihedral deconstructor
   */
  ~CrossDihedral(){}
  /**
   * Member operator = copies one CrossDihedral to the other
   */
  CrossDihedral &operator=(const CrossDihedral &);
  /**
   * Method to set the type of the CrossDihedral
   */
  void setType(int i){d_type = i;}
  /**
   * Accessor, returns the i-th atom in the cross dihedral as a const (0<i<3)
   */
  int operator[](int i)const{return d_a[i];}
  /**
   * Accessor, returns the i-th atom in the cross dihedral (0<i<3)
   */
  int &operator[](int i){return d_a[i];}
  /**
   * Accessor, returns the type of the CrossDihedral
   */
  int type()const{return d_type;}
};
/**
 * @relates CrossDihedral
 * Operator < compares two cross dihedral angles a and b
 *
 * a is defined to be smaller than b if <ol>
 * <li> the first dihedral a is smaller than the first dihedral of b
 * <li> or the second dihedral a is smaller than second dihedral of b
 *</ol>
 * @param a,b the CrossDihedral to be compared
 * @return 1 if a<b 0 otherwise
 */
int operator<(const CrossDihedral &a, const CrossDihedral &b);

}
#endif
