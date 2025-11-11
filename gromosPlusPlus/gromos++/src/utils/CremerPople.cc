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
#include "CremerPople.h"

#include <cassert>
#include <vector>
#include <cmath>
#include <cstdio>

#include "../gmath/Vec.h"
#include "../gmath/Physics.h"

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE leus

void util::cremerpople::calcCOM(const std::vector<gmath::Vec> &R, gmath::Vec &Rcm) {
  for (unsigned int i = 0; i < 3; ++i) {
    Rcm[i] = 0.0;
  }
  // vectors R' and R'' definitions
  for (unsigned int i = 0; i < 6; ++i) {
    Rcm += R[i];
  }
  Rcm /= 6.0;
}

void util::cremerpople::shiftToCOM(std::vector<gmath::Vec> &R, const gmath::Vec &Rcm) {
  for (unsigned int i = 0; i < 6; ++i) {
    R[i] -= Rcm;
  }
}

void util::cremerpople::calcPlane(const std::vector<gmath::Vec> &R, gmath::Vec &Rp, gmath::Vec &Rdp, gmath::Vec &RpxRdp) {
  for (unsigned int i = 0; i < 3; ++i) {
    Rp[i] = 0.0;
    Rdp[i] = 0.0;
  }
  // vectors R' and R'' definitions
  for (unsigned int i = 0; i < 6; ++i) {
    Rp += util::cremerpople::puckW[i] * R[i];
    Rdp += util::cremerpople::puckV[i] * R[i];
  }

  //normal vector of plane
  RpxRdp = Rp.cross(Rdp);
}

double util::cremerpople::calcZeta(const std::vector<gmath::Vec> &R, const gmath::Vec &RpxRdp, std::vector<double> &zeta) {
  double Q = 0.0;
  gmath::Vec n = RpxRdp / RpxRdp.abs();
  for (unsigned int i = 0; i < 6; ++i) {
    zeta[i] = n.dot(R[i]);
    Q += zeta[i] * zeta[i];
  }
  return sqrt(Q);
}

void util::cremerpople::calcGradZeta(const std::vector<gmath::Vec> &R, const gmath::Vec &Rp, const gmath::Vec &Rdp, const gmath::Vec &RpxRdp, const std::vector<double> &zeta, gmath::Vec &gizj, unsigned int index_i, unsigned int index_j) {
  gmath::Vec sum1, sum2, RdpxRj, RjxRp;
  double Ei = 0.0,
          Fi = 0.0;
  double RpRp, RpRdp, RdpRdp;
  double norm = RpxRdp.abs();
  for (unsigned int i = 0; i < 6; ++i) {
    Ei += Delta(index_i, i) * puckW[i];
    Fi += Delta(index_i, i) * puckV[i];
  }
  
  RpRp = Rp.dot(Rp);
  RdpRdp = Rdp.dot(Rdp);
  RpRdp = Rp.dot(Rdp);
  
  RdpxRj = Rdp.cross(R[index_j]);
  RjxRp = R[index_j].cross(Rp);
  sum1 = 2.0 * Ei * (Rp * RdpRdp - Rdp * RpRdp) + 2.0 * Fi * (Rdp * RpRp - Rp * RpRdp);
  sum2 = Delta(index_i, index_j) * RpxRdp + Ei * RdpxRj + Fi*RjxRp;
  
  gizj = -zeta[index_j] * sum1 / (2.0 * norm * norm) + sum2 / norm;
}

void util::cremerpople::calcGradQ(const std::vector<gmath::Vec> &R, const gmath::Vec &Rp, const gmath::Vec &Rdp, const gmath::Vec &RpxRdp, const std::vector<double> &zeta, gmath::Vec &giQ, unsigned int index) {
  double Q = 0.0;
  gmath::Vec gizj;
  for (unsigned int i = 0; i < 3; i++) {
    giQ[i] = 0.0;
  }

  for (unsigned int j = 0; j < 6; j++) {
    Q += zeta[j]*zeta[j];
    util::cremerpople::calcGradZeta(R, Rp, Rdp, RpxRdp, zeta, gizj, index, j);
    giQ += zeta[j] * gizj;
  }
  Q = sqrt(Q);
  giQ /= Q;
}

double util::cremerpople::calcPhi(const std::vector<double> &zeta) {
  double a = 0.0,
          b = 0.0;
  for (unsigned int i = 0; i < 6; ++i) {
    a += util::cremerpople::puckBIGW[i] * zeta[i];
    b += util::cremerpople::puckBIGV[i] * zeta[i];
  }
  
  double phi = atan(-a / b);
  if (b >= 0.0) {
    if (phi < 0.0) {
      phi += 2.0 * gmath::PhysConst().get_pi();
    }
  } else {
    phi += gmath::PhysConst().get_pi();
  }
  return phi;
}

void util::cremerpople::calcGradPhi(const std::vector<gmath::Vec> &R, const gmath::Vec &Rp, const gmath::Vec &Rdp, const gmath::Vec &RpxRdp, const std::vector<double> &zeta, gmath::Vec &giphi, unsigned int index) {
  gmath::Vec gizj, ga, gb;
  ga = 0.0;
  gb = 0.0;
  double a = 0.0,
          b = 0.0;
  for (unsigned int i = 0; i < 6; ++i) {
    util::cremerpople::calcGradZeta(R, Rp, Rdp, RpxRdp, zeta, gizj, index, i);
    a += util::cremerpople::puckBIGW[i] * zeta[i];
    b += util::cremerpople::puckBIGV[i] * zeta[i];
    ga += util::cremerpople::puckBIGW[i] * gizj;
    gb += util::cremerpople::puckBIGV[i] * gizj;
  }
  giphi = (a * gb - b * ga) / (a * a + b * b);
}

double util::cremerpople::calcTheta(const std::vector<double> &zeta) {
  double a = 0.0,
          b = 0.0,
          c = 0.0;
  for (unsigned int i = 0; i < 6; ++i) {
    a += util::cremerpople::puckBIGW[i] * zeta[i];
    b += util::cremerpople::puckBIGV[i] * zeta[i];
    c += util::cremerpople::puckSIGN[i] * zeta[i];
  }
  double theta = acos(c / sqrt(2. * (a * a + b * b) + c * c));
  return theta;
}

void util::cremerpople::calcGradTheta(const std::vector<gmath::Vec> &R, const gmath::Vec &Rp, const gmath::Vec &Rdp, const gmath::Vec &RpxRdp, const std::vector<double> &zeta, gmath::Vec &githeta, unsigned int index) {
  gmath::Vec gizj, ga, gb, gc;
  ga = 0.0;
  gb = 0.0;
  gc = 0.0;
  double a = 0.0,
          b = 0.0,
          c = 0.0,
          Q2 = 0.0;
  for (unsigned int i = 0; i < 6; ++i) {
    util::cremerpople::calcGradZeta(R, Rp, Rdp, RpxRdp, zeta, gizj, index, i);
    a += util::cremerpople::puckBIGW[i] * zeta[i];
    b += util::cremerpople::puckBIGV[i] * zeta[i];
    c += util::cremerpople::puckSIGN[i] * zeta[i];
    ga += util::cremerpople::puckBIGW[i] * gizj;
    gb += util::cremerpople::puckBIGV[i] * gizj;
    gc += util::cremerpople::puckSIGN[i] * gizj;
    Q2 += zeta[i] * zeta[i];
  }
    
  double sqrta2b2 = sqrt(a * a + b * b);
  githeta = c * (a * ga + b * gb) / sqrta2b2 - sqrta2b2 * gc;
  githeta = githeta / (3. * sqrt(2.) * Q2);
}
