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

// gcore_Box.h

#ifndef BOX_H
#define BOX_H


#include <cassert>
#include <vector>

#include "../gmath/Vec.h"
#include "../gromos/Exception.h"


namespace gmath {
  class Vec;
}

namespace gcore {

  /**
   * Class Box
   * Purpose: contains the box size and angle
   *
   * Description:
   * The Box class contains the box dimensions and angle
   *
   * @class Box
   * @author R. Buergi
   * @ingroup gcore
   * @sa gcore::System
   */
  class Box {
  public:

    enum boxshape_enum {
      vacuum = 0, rectangular = 1, triclinic = 2, truncoct = -1
    };

    enum boxformat_enum {
      box96, triclinicbox, genbox
    };

  private:
    std::vector<gmath::Vec> d_dim;
    double d_K_L_M;
    gmath::Vec d_origin;
    std::vector<gmath::Vec> d_cross_K_L_M;
    boxshape_enum d_ntb;
    boxformat_enum d_boxformat;

  public:
    // constructurs

    /**
     * Box constructor
     * @param x, y, z box dimensions
     */
    Box(double x = 0, double y = 0, double z = 0) :
    d_dim(3), d_K_L_M(0), d_origin(0.0),
    d_cross_K_L_M(3), d_ntb(vacuum), d_boxformat(box96) {
      d_dim[0][0] = x;
      d_dim[1][1] = y;
      d_dim[2][2] = z;
    }

    /**
     * Box constructor
     * @param K, L, M box vectors
     */
    Box(gmath::Vec K, gmath::Vec L, gmath::Vec M) :
    d_dim(3), d_K_L_M(0), d_origin(0.0),
    d_cross_K_L_M(3), d_ntb(vacuum), d_boxformat(triclinicbox) {
      d_dim[0] = K;
      d_dim[1] = L;
      d_dim[2] = M;
    }

    /**
     * Box constructor from dimensions and Euler angles and box rotation
     * @param bound boundary conditions
     * @param a, b, c box dimensions
     * @param alpha, beta, gamma Euler angles
     * @param phi, theta, psi box rotation angles
     */
    Box(boxshape_enum bound,
            double a, double b, double c,
            double alpha, double beta, double gamma,
            double phi, double theta, double psi, double X = 0.0, double Y = 0.0, double Z = 0.0);

    /**
     * Box copy constructor
     * @param b Box to be copied
     */
    Box(const Box&b) :
    d_dim(b.d_dim), d_K_L_M(b.d_K_L_M), d_origin(0.0),
    d_cross_K_L_M(b.d_cross_K_L_M), d_ntb(b.d_ntb),
    d_boxformat(b.d_boxformat) {
    }
    /**
     * Update volume and cross product for the triclinic box
     */
    void update_triclinic();

    // accessors
    Box & operator=(const Box&);

    /**
     * Accessor, return the K-vector of a generalized box
     */
    gmath::Vec &K();
    /**
     * Accessor, return the K-vector of a generalized box
     */
    gmath::Vec const & K()const;
    /**
     * Accessor, return the L-vector of a generalized box
     */
    gmath::Vec &L();
    /**
     * Accessor, return the L-vector of a generalized box
     */
    gmath::Vec const & L()const;
    /**
     * Accessor, return the M-vector of a generalized box
     */
    gmath::Vec &M();
    /**
     * Accessor, return the M-vector of a generalized box
     */
    gmath::Vec const & M()const;
    /**
     * Sets the length of the K vector
     */
    void stretch_K(double l);
    /**
     * Sets the length of the L vector
     */
    void stretch_L(double l);
    /**
     * Sets the length of the M vector
     */
    void stretch_M(double l);
    /**
     * Accessor, return the value of ntb
     */
    void setNtb(boxshape_enum b);
    /**
     * Accessor, return the vlaue of ntb
     */
    boxshape_enum ntb()const;
    /**
     * Accessor, return the format of the box
     */
    boxformat_enum &boxformat();
    /**
     * Accessor, return the format of the box
     */
    boxformat_enum boxformat()const;

    /**
     * Accessor, return the (triclinic) volume of a generalized box
     */
    double K_L_M();
    /**
     * Accessor, return the (triclinic) volume of a generalized box
     */
    double K_L_M()const;
    /**
     * Accessor, return the weird cross product of a generalized box
     */
    std::vector<gmath::Vec> & cross_K_L_M();
    /**
     * Accessor, return the weird cross product of a generalized box
     */
    std::vector<gmath::Vec> const & cross_K_L_M()const;
    /**
     * returns the alpha angle in degree
     */
    double alpha()const;
    /**
     * returns the beta angle in degree
     */
    double beta()const;
    /**
     * returns the gamma angle in degree
     */
    double gamma()const;
    /**
     * returns the X value
     */
    double X()const;
    /**
     * returns the Y value
     */
    double Y()const;
    /**
     * returns the Z value
     */
    double Z()const;
  };

  inline gcore::Box &Box::operator=(gcore::Box const &b) {
    if (this != &b) {
      d_dim = b.d_dim;
      d_K_L_M = b.d_K_L_M;
      d_cross_K_L_M = b.d_cross_K_L_M;
      d_ntb = b.d_ntb;
      d_boxformat = b.d_boxformat;
      d_origin = b.d_origin;
    }
    return *this;
  }

  inline gmath::Vec &Box::K() {
    return d_dim[0];
  }

  inline gmath::Vec const & Box::K()const {
    return d_dim[0];
  }

  inline gmath::Vec &Box::L() {
    return d_dim[1];
  }

  inline gmath::Vec const & Box::L()const {
    return d_dim[1];
  }

  inline gmath::Vec &Box::M() {
    return d_dim[2];
  }

  inline gmath::Vec const & Box::M()const {
    return d_dim[2];
  }

  inline void Box::stretch_K(double l) {
    double l_ = d_dim[0].abs();
    if (l_ == 0)
      gromos::Exception("Box.h", "Cannot stretch a vector of length 0!");
    d_dim[0] *= l / l_;
  }

  inline void Box::stretch_L(double l) {
    double l_ = d_dim[1].abs();
    if (l_ == 0)
      gromos::Exception("Box.h", "Cannot stretch a vector of length 0!");
    d_dim[1] *= l / l_;
  }

  inline void Box::stretch_M(double l) {
    double l_ = d_dim[2].abs();
    if (l_ == 0)
      gromos::Exception("Box.h", "Cannot stretch a vector of length 0!");
    d_dim[2] *= l / l_;
  }

  inline Box::boxshape_enum Box::ntb()const {
    return d_ntb;
  }

  inline Box::boxformat_enum & Box::boxformat() {
    return d_boxformat;
  }

  inline Box::boxformat_enum Box::boxformat()const {
    return d_boxformat;
  }

  inline double Box::K_L_M() {
    return d_K_L_M;
  }

  inline double Box::K_L_M()const {
    return d_K_L_M;
  }

  inline std::vector<gmath::Vec> & Box::cross_K_L_M() {
    return d_cross_K_L_M;
  }

  inline std::vector<gmath::Vec> const & Box::cross_K_L_M()const {
    return d_cross_K_L_M;
  }

  inline double Box::X()const {
    return d_origin[0];
  }

  inline double Box::Y()const {
    return d_origin[1];
  }

  inline double Box::Z()const {
    return d_origin[2];
  }

} /*namespace*/
#endif
