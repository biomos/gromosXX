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

// VirtualAtomType.h
#ifndef INCLUDED_VIRTUALATOMTYPE
#define INCLUDED_VIRTUALATOMTYPE

namespace gcore{

  /**
   * Class VirtualAtomType
   * Purpose: contains a virtual atom type
   *
   * Description:
   * Contains two distances that may be relevant to calculate the position of a virtual atom
   *
   * @class VirtualAtomType
   * @author R. Buergi, C. Oostenbrink
   * @ingroup gcore
   * @sa gcore::GromosForceField
   */

class VirtualAtomType
{
  int d_code;
  double d_dis1;
  double d_dis2;
 public:
  /**
   * VirtualAtomType constructor
   * @param c    bond code
   * @param dis1 distance 1 (DISH)
   * @param dis2 distance 2 (DISC)
   */
  VirtualAtomType(int c=0, double dis1=0, double dis2=0): d_code(c), d_dis1(dis1), d_dis2(dis2){}
  /**
   * VirtualAtomType copyconstructor
   * @param va VirtualAtomType to be copied
   */
  VirtualAtomType(const VirtualAtomType& va): d_code(va.d_code), d_dis1(va.d_dis1), d_dis2(va.d_dis2){}
  /** 
   * Member operator=, assign distances of one virtual atom type to another
   */
  VirtualAtomType &operator=(const VirtualAtomType &va);
  /**
   * BondType deconstuctor
   */
  ~VirtualAtomType(){}
  /**
   * Accessor, returns the integer code
   */
  int code()const{return d_code;}
  /**
   * Accessor, returns the first distance
   */
  double dis1()const{return d_dis1;}
  /**
   * Accessor, returns the second distance
   */
  double dis2()const{return d_dis2;}
};

}
#endif



