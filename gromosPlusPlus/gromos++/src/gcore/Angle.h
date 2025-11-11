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

// gcore_Angle.h

#ifndef INCLUDED_GCORE_ANGLE
#define INCLUDED_GCORE_ANGLE

namespace gcore{

  /**
   * Class Angle
   * Purpose: contains a gromos96 angle
   *
   * Description:
   * Contains the atoms and type making up an angle. The atoms are sorted
   * in such a way that the atom number of the first atom is lower than the
   * atom number of the last atom.
   *
   * @class Angle
   * @version $Date: Dec 11 2001
   * @author R. Buergi
   * @ingroup gcore
   * @sa gcore::AngleType
   * @sa gcore::MoleculeTopology
   */
class Angle{
  int d_a[3];
  int d_type;
  // not implemented
  Angle();
 public:
  /**
   * Angle constructor
   * Constructs a new angle defined by atoms a, b and c. The angle is
   * stored in such a way that a < c. 
   * @param a,b,c atomnumbers of the atoms involved
   */
  Angle(int a, int b, int c, bool warn = true);
  /**
   * Angle copy-constructor
   * Construct a new angle and copy the specified angle
   * @param & angle to be copied
   */
  Angle(const Angle &);
  /**
   * Angle deconstructor
   */
  ~Angle(){}
  /**
   * Method setType, sets the angle type of the angle
   * @param i angle type to be set
   */
  void setType(int i){d_type=i;}
  /**
   * Member operator =, assigns atoms and angletype of one angle
   * to the next.
   */
  Angle &operator=(const Angle &);
  /**
   * Member operator [], returns the i-th atom forming the angle as a const
   * @param i a number from 0 to 2
   * @return The atom number of the i-th atom in the angle
   */
  int operator[](int i)const{return d_a[i];}
  /**
   * Member operator [], returns the i-th atom forming the angle
   * @param i a number from 0 to 2
   * @return The atom number of the i-th atom in the angle
   */
  int &operator[](int i){return d_a[i];}
  /**
   * Accessor type, returns the angle type of the angle
   * @return The angle type of the angle
   */
  int type()const{return d_type;}
};

/**
 * @relates Angle
 * Operator < compares two angles a and b
 *
 * a is defined to be smaller than b if<br><ol>
 * <li> Atomnumber of the centre atom of a is lower than the atomnumber of 
 *      the central atom of b (a[1] < b[1]), or
 * <li> The centre atoms have the same atomnumber, but the atomnumber of 
 *      the first atom in a is smaller than that in b (a[1]==b[1] && 
 *      a[0]<b[0]), or
 * <li> The first two atoms in a and b are the same, but the third atom in
 *      a has a lower number than that in b (a[0]==b[0] && a[1]==b[1] && 
 *      a[2]<b[2]), or
 * <li> All atoms are the same, but a has a lower angle type than b 
 *      (a[0]==b[0] && a[1]==b[1] && a[2]==b[2] && a.type() < b.type())
 * </ol>
 * @param a, b angles to be compared
 * @return 0 if a > b and 1 if a<b
 */
int operator<(const Angle &a, const Angle &b);

}
#endif

