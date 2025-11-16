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

// AngleType.h
#ifndef INCLUDED_ANGLETYPE
#define INCLUDED_ANGLETYPE

namespace gcore{

  /**
   * Class AngleType
   * Purpose: contains a covalent bond angle type
   *
   * Description:
   * Contains the optimum angle and force constant for a covalent bond
   * angle. The potential energy for an angle is defined as
   * @f[ V^{angle}=\frac{1}{2}Kc_{\theta_n}\left[\cos{\theta_n} - \cos{\theta_{0_n}}\right]^2@f]
   * with @f$Kc_{\theta_n}@f$ in kJ/mol, or as
   * @f[ V^{angle}=\frac{1}{2}Ka_{\theta_n}\left[\theta_n - \theta_{0_n}\right]^2@f]
   * with @f$Ka_{\theta_n}@f$ in kJ/mol/deg^2
   *
   * @class AngleType
   * @author R. Buergi, D. Geerke
   * @ingroup gcore
   * @sa gcore::Angle
   * @sa gcore::GromosForceField
   */
class AngleType
{
  int d_code;
  double d_t0;
  double d_fc;
  double d_afc;
 public:
  /**
   * AngleType constructor
   * @param c   integer code of the angle type
   * @param fc  force constant (@f$Kc_{\theta_n}@f$)
   * @param afc force constant (@f$Ka_{\theta_n}@f$)
   * @param l   optimum angle (@f$\theta_{0_n}@f$)
   */
  AngleType(int c, double fc, double afc, double l): d_code(c), d_t0(l), 
          d_fc(fc), d_afc(afc) {}
  /**
   * AngleType constructor
   * The hamonic force constant @f$Ka_{\theta_n}@f$ is calculated as described in
   * volume 2.
   * @param c   integer code of the angle type
   * @param fc  force constant (@f$Kc_{\theta_n}@f$)
   * @param l   optimum angle (@f$\theta_{0_n}@f$)
   */
  AngleType(int c=0, double fc=0, double l=0);
  /**
   * AngleType copyconstructor
   * @param b AngleType to be copied
   */
  AngleType(const AngleType& b):d_code(b.d_code), d_t0(b.d_t0), d_fc(b.d_fc),
          d_afc(b.d_afc){}
  /** 
   * Member operator=, assign force constant and optimum angle of one
   * AngleType to the other
   */
  AngleType &operator=(const AngleType &b);
  /**
   * AngleType deconstuctor
   */
  ~AngleType(){}
  /**
   * Accessor, returns the integer code
   */
  int code()const{return d_code;}
  /**
   * Accessor, returns the optimum angle (@f$\theta_{0_n}@f$)
   */
  double t0()const{return d_t0;}
  /**
   * Accessor, returns the force constant (@f$Kc_{\theta_n}@f$)
   */
  double fc()const{return d_fc;}
  /**
   * Accessor, returns the force constant (@f$Ka_{\theta_n}@f$)
   */
  double afc()const{return d_afc;}
};

}
#endif
