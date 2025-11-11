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

// gcore_Box.cc

#include "Box.h"

#include <cassert>
#include <cmath>
#include <math.h>

#include "../gromos/Exception.h"
#include "../gmath/Vec.h"
#include "../gmath/Matrix.h"

void gcore::Box::update_triclinic() {
  d_K_L_M = (K().cross(L())).dot(M());
  if (d_K_L_M == 0.0) return; // weird box, vacuum or similar
  d_cross_K_L_M[0] = L().cross(M()) / -d_K_L_M;
  d_cross_K_L_M[1] = K().cross(M()) / d_K_L_M;
  d_cross_K_L_M[2] = K().cross(L()) / -d_K_L_M;
}

void gcore::Box::setNtb(boxshape_enum b) {
  d_ntb = b;
}

gcore::Box::Box(gcore::Box::boxshape_enum bound, double a, double b, double c, double alpha, double beta, double gamma,
        double phi, double theta, double psi, double X, double Y, double Z) {
  d_boxformat = gcore::Box::genbox;

  d_dim.resize(3);
  d_cross_K_L_M.resize(3);
  d_ntb = bound;

  d_origin = gmath::Vec(X,Y,Z);

  if (d_ntb == gcore::Box::vacuum) {

    K() = gmath::Vec(0.0, 0.0, 0.0);
    L() = gmath::Vec(0.0, 0.0, 0.0);
    M() = gmath::Vec(0.0, 0.0, 0.0);
    update_triclinic();

    return;
  }

  if (d_ntb == gcore::Box::rectangular || d_ntb == gcore::Box::truncoct) {
    if (alpha != 90.0 || beta != 90.0 || gamma != 90.0)
      throw gromos::Exception("GENBOX", "For rectangular and truncated octahedral boxes, alpha, beta"
            " and gamma should be 90 degrees");
    if (phi != 0.0 || theta != 0.0 || psi != 0.0)
      throw gromos::Exception("GENBOX","For rectangular and truncated octahedral boxes, phi, theta"
            " and phi should be 0 degrees");
    K() = gmath::Vec(a, 0.0, 0.0);
    L() = gmath::Vec(0.0, b, 0.0);
    M() = gmath::Vec(0.0, 0.0, c);
    update_triclinic();
    return;
  }

  // lets generate the L_ = R_ * S_ (L_ = (K, L, M)) according to
  // Christen et al. (2005): J.Comp.Chem., Vol. 26, No. 16, 1719 - 1751.
  alpha *= M_PI / 180.0;
  beta *= M_PI / 180.0;
  gamma *= M_PI / 180.0;
  const double cosalpha = cos(alpha);
  const double cosbeta = cos(beta);
  const double cosgamma = cos(gamma);
  const double sinbeta = sin(beta);
  const double singamma = sin(gamma);
  const double cosdelta = (cosalpha - cosbeta * cosgamma) / (sinbeta * singamma);
  const double sindelta = sqrt(1 - cosdelta * cosdelta);
  // the three (columns) vectors of the transformation matrix S_ = (S_1, S_2, S_3)
  gmath::Vec S_1(a, 0.0, 0.0);
  gmath::Vec S_2(b * cosgamma, b * singamma, 0.0);
  gmath::Vec S_3(c * cosbeta, c * sinbeta * cosdelta, c * sinbeta * sindelta);
  gmath::Matrix S_(S_1, S_2, S_3);

  // now the rotation matrix R_
  // if no rotation has to be applied: set K = S_1, L = S_2 and M = S3
  if (phi != 0 || theta != 0 || psi != 0) {
    phi *= M_PI / 180.0;
    theta *= M_PI / 180.0;
    psi *= M_PI / 180.0;
    const double cosphi = cos(phi);
    const double costheta = cos(theta);
    const double cospsi = cos(psi);
    const double sinphi = sin(phi);
    const double sintheta = sin(theta);
    const double sinpsi = sin(psi);
    // the three (columns) vectors of the transformation matrix R_ = (R_1, R_2, R_3)
    gmath::Vec R_1(costheta * cosphi, costheta * sinphi, -sintheta);
    gmath::Vec R_2(sinpsi * sintheta * cosphi - cospsi * sinphi,
            sinpsi * sintheta * sinphi + cospsi * cosphi,
            sinpsi * costheta);
    gmath::Vec R_3(cospsi * sintheta * cosphi + sinpsi * sinphi,
            cospsi * sintheta * sinphi - sinpsi * cosphi,
            cospsi * costheta);
    gmath::Matrix R_(R_1, R_2, R_3);
    gmath::Matrix L_(R_*S_);
    K() = gmath::Vec(L_(0,0), L_(1,0), L_(2,0));
    L() = gmath::Vec(L_(0,1), L_(1,1), L_(2,1));
    M() = gmath::Vec(L_(0,2), L_(1,2), L_(2,2));
    update_triclinic();
  } else {
    K() = S_1;
    L() = S_2;
    M() = S_3;
    update_triclinic();
    return;
  }
  return;
}

double gcore::Box::alpha() const {
  return acos(L().dot(M())/(L().abs()*M().abs()))*180.0/M_PI;
}

double gcore::Box::beta() const {
  return acos(K().dot(M())/(K().abs()*M().abs()))*180.0/M_PI;
}

double gcore::Box::gamma() const {
  return acos(K().dot(L())/(K().abs()*L().abs()))*180.0/M_PI;
}

