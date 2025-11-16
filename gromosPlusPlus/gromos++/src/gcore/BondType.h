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

// BondType.h
#ifndef INCLUDED_BONDTYPE
#define INCLUDED_BONDTYPE

namespace gcore{

  /**
   * Class BondType
   * Purpose: contains a covalent bond type
   *
   * Description:
   * Contains the optimum bond length (@f$b_{0_n}@f$), quartic force constant
   * (@f$Kq_{b_n}@f$) and a harmonic force constant (@f$Kh_{b_n}@f$) for a covalent bond.
   * The potential energy for a bond is defined as
   * @f[ V^{bond}=\frac{1}{4}Kh_{b_n}\left[b_n^2 - b_{0_n}^2\right]^2@f]
   * or
   * @f[ V^{bond}=\frac{1}{2}Kq_{b_n}\left[b_n^1 - b_{0_n}^1\right]^2@f]
   *
   * @class BondType
   * @author R. Buergi, D. Geerke
   * @ingroup gcore
   * @sa gcore::Bond
   * @sa gcore::GromosForceField
   */

class BondType
{
  int d_code;
  double d_fc;
  double d_hfc;
  double d_b0;
 public:
  /**
   * BondType constructor
   * @param c   bond code
   * @param fc  quartic force constant  (@f$Kq_{b_n}@f$)
   * @param hfc harmonic force constant (@f$Kh_{b_n}@f$)
   * @param l   equilibrium bond length (@f$b_{0_n}@f$)
   */
  BondType(int c, double fc, double hfc, double l): d_code(c), d_fc(fc), d_hfc(hfc), d_b0(l){}
  /**
   * BondType constructor
   * The harmonic force constant @f$Kh_{b_n}@f$ is calculated
   * @f$[Kh_{b_n} = 2b_{0_n}^2 Kq_{b_n}]@f$
   * @param c   bond code
   * @param fc quartic force constant (@f$Kq_{b_n}@f$)
   * @param l   equilibrium bond length (@f$b_{0_n}@f$)
   */
  BondType(int c=0, double fc=0, double l=0, bool quartic=true);
  //BondType(int c=0, double fc=0, double l=0);
  /**
   * BondType copyconstructor
   * @param b BondType to be copied
   */
  BondType(const BondType& b): d_code(b.d_code), d_fc(b.d_fc), d_hfc(b.d_hfc),
          d_b0(b.d_b0){}
  /** 
   * Member operator=, assign force constant and optimum bond lenght of one
   * BondType to the other
   */
  BondType &operator=(const BondType &b);
  /**
   * BondType deconstuctor
   */
  ~BondType(){}
  /**
   * Accessor, returns the integer code
   */
  int code()const{return d_code;}
  /**
   * Accessor, returns the optimum bond length (@f$b_{0_n}@f$)
   */
  double b0()const{return d_b0;}
  /**
   * Accessor, returns the quartic force constant (@f$Kq_{b_n}@f$)
   */
  double fc()const{return d_fc;}
  /**
   * Accessor, returns the harmonic force constant (@f$Kh_{b_n}@f$)
   */
  double hfc()const{return d_hfc;}
};

}
#endif



