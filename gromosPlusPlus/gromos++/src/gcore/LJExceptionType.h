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

// LJExceptionType.h
#ifndef INCLUDED_LJEXCEPTIONTYPE
#define INCLUDED_LJEXCEPTIONTYPE

namespace gcore{

  /**
   * Class LJExceptionType
   * Purpose: contains a LJException type
   *
   * Description:
   * Contains the C12 and C6 parameters of a LJ exception type as well as its
   * type number
   *
   * @class LJExceptionType
   * @author A.P. Eichenberger
   * @ingroup gcore
   */

class LJExceptionType
{
  int d_code;
  double d_c12;
  double d_c6;
 public:
  /**
   * LJExceptionType constructor
   */
  LJExceptionType(int c=-1, double c12=0.0, double c6=0.0):d_code(c), d_c12(c12), d_c6(c6) {};
  /**
   * LJExceptionType copyconstructor
   * @param b LJExceptionType to be copied
   */
  LJExceptionType(const LJExceptionType& b): d_code(b.d_code), d_c12(b.d_c12), d_c6(b.d_c6){}
  /** 
   * Member operator =
   */
  LJExceptionType &operator=(const LJExceptionType &b);
  /**
   * LJExceptionType deconstuctor
   */
  ~LJExceptionType(){}
  /**
   * Accessor, returns the integer code
   */
  int code()const{return d_code;}
  /**
   * Accessor, returns the C12 parameter
   */
  double c12()const{return d_c12;}
  /**
   * Accessor, returns the C6 parameter
   */
  double c6()const{return d_c6;}
};

}
#endif



