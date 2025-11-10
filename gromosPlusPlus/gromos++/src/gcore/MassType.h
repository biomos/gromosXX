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

// MassType.h
#ifndef INCLUDED_MASSTYPE
#define INCLUDED_MASSTYPE

namespace gcore{
  /**
   * Class MassType
   * Purpose: contains a gromos96 mass type
   *
   * Description:
   * The gromos96 integer mass types are contained in the class...
   * MassType!
   *
   * @class MassType
   * @author C. Oostenbrink
   * @ingroup gcore
   * @sa gcore::GromosForceField
   */
class MassType
{
  int d_n;
  double d_am;
 public:
  /**
   * MassType constructor
   * @param n The gromos96 integer mass code
   * @param l The corresponding mass in atomic units
   */
  MassType(int n=0, double l=0): d_n(n), d_am(l){}
  /**
   * MassType copy constructor
   * @param b The MassType to be copied
   */
  MassType(const MassType& b):d_n(b.d_n), d_am(b.d_am){}
  /**
   * Member operator = copies one MassType into the other
   */
  MassType &operator=(const MassType &b);
  /**
   * MassType deconstructor
   */
  ~MassType(){}
  /**
   * Accessor, returns the Integer Mass Code
   */
  int n()const{return d_n;}
  /**
   * Accessor, returns the atomic mass in atomic units
   */
  double am()const{return d_am;}
};

}
#endif

