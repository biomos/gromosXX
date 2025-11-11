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
#include "Gch.h"

#include <cmath>
#include <cstdlib>
#include <math.h>
#include <sstream>
#include <vector>
#include <string>
#include <cassert>
#include <iostream>

#include "../gio/InTopology.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/Bond.h"
#include "../gcore/BondType.h"
#include "../gcore/Angle.h"
#include "../gcore/AngleType.h"
#include "../gcore/Dihedral.h"
#include "../gromos/Exception.h"
#include "../gmath/Vec.h"
#include "../gmath/Matrix.h"
#include "../utils/Neighbours.h"
#include "../gcore/GromosForceField.h"

using namespace gcore;
using namespace std;

double find_bond(gcore::System &sys, gcore::GromosForceField &gff, int m, gcore::Bond b, double guess);
double find_angle(gcore::System &sys, gcore::GromosForceField &gff, int m, gcore::Angle a, double guess);
int find_dihedral(gcore::System &sys, int m, int i, int j, std::vector<int> h);

int utils::generate_hcoordinates(System &sys, GromosForceField &gff, int m, int a, std::vector<int> &h, std::vector<int> &nh, int geom, double eps) {
  int count = 0;
  switch (geom) {
    case(1):
    {

      double bond = find_bond(sys, gff, m, Bond(a, h[0], false), 0.1);
      Vec v0 = sys.mol(m).pos(a) - sys.mol(m).pos(h[0]);
      if (fabs(v0.abs() - bond) / bond > eps) {

        double angle = M_PI / 180.0 * find_angle(sys, gff, m,
                Angle(nh[0], a, h[0]), 109.5);

        // in case of ARGN, HH1 will be in cis-conformation instead of anti
        // (dihedral NE-CZ-NH1-HH1)
        // Claudio M. Soares, private communication, May 2011
        if (sys.mol(m).topology().resName(sys.mol(m).topology().resNum(nh[0])) == "ARGN") {
          angle *= -1;
        }

        int fourth = find_dihedral(sys, m, nh[0], a, h);

        Vec v1 = (sys.mol(m).pos(fourth) - sys.mol(m).pos(nh[0])).normalize();
        Vec v2 = (sys.mol(m).pos(nh[0]) - sys.mol(m).pos(a)).normalize();
        Vec v4 = (v1.cross(v2)).normalize();
        Vec v5 = (v2.cross(v4)).normalize();
        Vec v6 = bond * cos(angle) * v2 - bond * sin(angle) * v5;

        sys.mol(m).pos(h[0]) = sys.mol(m).pos(a) + v6;
        count++;
      }
      break;
    }

    case(2):
    {
      double bond = find_bond(sys, gff, m, Bond(a, h[0], false), 0.1);
      Vec v0 = sys.mol(m).pos(a) - sys.mol(m).pos(h[0]);
      if (fabs(v0.abs() - bond) / bond > eps) {

        Vec v1 = (sys.mol(m).pos(nh[0]) - sys.mol(m).pos(a)).normalize();
        Vec v2 = (sys.mol(m).pos(nh[1]) - sys.mol(m).pos(a)).normalize();
        Vec v3 = (v1 + v2).normalize();
        sys.mol(m).pos(h[0]) = sys.mol(m).pos(a) - bond*v3;
        count++;
      }
      break;
    }
    case(3):
    {
      double bond1 = find_bond(sys, gff, m, Bond(a, h[0], false), 0.1);
      double bond2 = find_bond(sys, gff, m, Bond(a, h[1], false), 0.1);
      Vec v01 = sys.mol(m).pos(a) - sys.mol(m).pos(h[0]);
      Vec v02 = sys.mol(m).pos(a) - sys.mol(m).pos(h[1]);
      if (fabs(v01.abs() - bond1) / bond1 > eps ||
              fabs(v02.abs() - bond2) / bond2 > eps) {

        double angle1 = M_PI / 180.0 * find_angle(sys, gff, m,
                Angle(nh[0], a, h[0], false), 120);
        double angle2 = M_PI / 180.0 * find_angle(sys, gff, m,
                Angle(nh[0], a, h[1], false), 120);
        double angle3 = M_PI / 180.0 * find_angle(sys, gff, m,
                Angle(h[0], a, h[1], false), 120);
        int fourth = find_dihedral(sys, m, nh[0], a, h);

        Vec v1 = (sys.mol(m).pos(fourth) - sys.mol(m).pos(nh[0])).normalize();
        Vec v2 = (sys.mol(m).pos(nh[0]) - sys.mol(m).pos(a)).normalize();
        Vec v4 = v1.cross(v2).normalize();
        Vec v5 = v2.cross(v4).normalize();
        Vec v6 = bond1 * cos(angle1) * v2
                - bond1 * sin(angle1) * v5;
        double A = bond2 * cos(angle2);
        double B = (bond1 * bond2 * cos(angle3) - A * v2.dot(v6)) / v5.dot(v6);
        double C = sqrt(bond2 * bond2 - A * A - B * B);
        Vec v7 = A * v2 + B * v5 + C * v4;

        // we do not want to consider if the individual hydrogens
        // need to be replaced, because if h[0] has a correct bond
        // length, but is placed at the position where we predict h[1]
        // this leads to crashes later on. If one was not good, we 
        // replace both
        sys.mol(m).pos(h[0]) = sys.mol(m).pos(a) + v6;
        count++;
        sys.mol(m).pos(h[1]) = sys.mol(m).pos(a) + v7;
        count++;
      }
      break;
    }
    case(4):
    {
      // very similar to the non-planar type of above
      double bond1 = find_bond(sys, gff, m, Bond(a, h[0], false), 0.1);
      double bond2 = find_bond(sys, gff, m, Bond(a, h[1], false), 0.1);
      double bond3 = find_bond(sys, gff, m, Bond(a, h[2], false), 0.1);
      Vec v01 = sys.mol(m).pos(a) - sys.mol(m).pos(h[0]);
      Vec v02 = sys.mol(m).pos(a) - sys.mol(m).pos(h[1]);
      Vec v03 = sys.mol(m).pos(a) - sys.mol(m).pos(h[2]);

      if (fabs(v01.abs() - bond1) / bond1 > eps ||
              fabs(v02.abs() - bond2) / bond2 > eps ||
              fabs(v03.abs() - bond3) / bond3 > eps) {

        double angle1 = M_PI / 180.0 * find_angle(sys, gff, m,
                Angle(nh[0], a, h[0], false), 109.5);
        double angle2 = M_PI / 180.0 * find_angle(sys, gff, m,
                Angle(nh[0], a, h[1], false), 109.5);
        double angle3 = M_PI / 180.0 * find_angle(sys, gff, m,
                Angle(nh[0], a, h[2], false), 109.5);
        double angle4 = M_PI / 180.0 * find_angle(sys, gff, m,
                Angle(h[0], a, h[1], false), 109.5);
        double angle5 = M_PI / 180.0 * find_angle(sys, gff, m,
                Angle(h[0], a, h[2], false), 109.5);
        int fourth = find_dihedral(sys, m, nh[0], a, h);

        Vec v1 = (sys.mol(m).pos(fourth) - sys.mol(m).pos(nh[0])).normalize();
        Vec v2 = (sys.mol(m).pos(nh[0]) - sys.mol(m).pos(a)).normalize();
        Vec v3 = (sys.mol(m).pos(a) - sys.mol(m).pos(h[0])).normalize();
        Vec v4 = v2.cross(v1);
        if(v4.abs() < 1e-10) {
          if(v2[0] != 0.0 || v2[1] != 0.0) {
            v4[0] = -v2[1];
            v4[1] = v2[0];
            v4[2] = 0;
          } else {
            v4[0] = 0;
            v4[1] = -v2[2];
            v4[2] = v2[1];
          }
          v4 = v4.normalize();
        } else {
          v4 = v4.normalize();
        }
        Vec v5 = v2.cross(v4).normalize();
        Vec v6 = bond1 * cos(angle1) * v2 - bond1 * sin(angle1) * v5;

        double A = bond2 * cos(angle2);
        double B = (bond1 * bond2 * cos(angle4) - A * v2.dot(v6)) / v5.dot(v6);
        double C = sqrt(bond2 * bond2 - A * A - B * B);
        Vec v7 = A * v2 + B * v5 + C * v4;
        A = bond3 * cos(angle3);
        B = (bond1 * bond3 * cos(angle5) - A * v2.dot(v6)) / v5.dot(v6);
        C = sqrt(bond3 * bond3 - A * A - B * B);
        Vec v8 = A * v2 + B * v5 - C * v4;

        // we do not want to consider if the individual hydrogens
        // need to be replaced, because if h[0] has a correct bond
        // length, but is placed at the position where we predict h[1]
        // this leads to crashes later on. If one was not good, we 
        // replace all three
        sys.mol(m).pos(h[0]) = sys.mol(m).pos(a) + v6;
        count++;
        sys.mol(m).pos(h[1]) = sys.mol(m).pos(a) + v7;
        count++;
        sys.mol(m).pos(h[2]) = sys.mol(m).pos(a) + v8;
        count++;
      }
      break;
    }
    case(5):
    {
      // likely to be a water molecule. Here we have to come up with some
      // random orientation.
      double bond1 = find_bond(sys, gff, m, Bond(a, h[0], false), 0.1);
      double bond2 = find_bond(sys, gff, m, Bond(a, h[1], false), 0.1);
      Vec v01 = sys.mol(m).pos(a) - sys.mol(m).pos(h[0]);
      Vec v02 = sys.mol(m).pos(a) - sys.mol(m).pos(h[1]);
      if (fabs(v01.abs() - bond1) / bond1 > eps ||
              fabs(v02.abs() - bond2) / bond2 > eps) {
        // first generate a standard molecule. If it really is a water 
        // molecule, there is probably no angle defined, but putting them 
        // at 109.5 degrees is not so bad.
        double angle = M_PI / 180.0 * find_angle(sys, gff, m,
                Angle(h[0], a, h[1], false), 109.5);
        Vec v1(0.0, 0.0, bond1);
        Vec v2(0.0, bond2 * sin(angle), bond2 * cos(angle));

        //get three random numbers for the angles
        //calculate sin and cos of these angles

        Vec angle_cos, angle_sin;
        for (int i = 0; i < 3; i++) {
          double ang = 2.0 * M_PI / RAND_MAX * double(rand());
          angle_cos[i] = cos(ang);
          angle_sin[i] = sin(ang);
        }

        // prepare a matrix to perform three random rotations  
        // The product of three rotation matrices about three axes
        /*
         * (  1.0   0.0   0.0)   ( cosy   0.0   siny)   ( cosx  -sinx   0.0)
         * (  0.0  cosz -sinz) X (  0.0   1.0    0.0) X ( sinx   cosx   0.0)
         * (  0.0  sinz  cosz)   (-siny   0.0   cosy)   (  0.0    0.0   1.0)
         */
        gmath::Matrix rot(3, 3);
        rot(0, 0) = angle_cos[0] * angle_cos[1];
        rot(1, 0) = angle_sin[0] * angle_cos[2]
                + angle_cos[0] * angle_sin[1] * angle_sin[2];
        rot(2, 0) = angle_sin[0] * angle_sin[2]
                - angle_cos[0] * angle_sin[1] * angle_cos[2];
        rot(0, 1) = -angle_sin[0] * angle_cos[1];
        rot(1, 1) = angle_cos[0] * angle_cos[2]
                - angle_sin[0] * angle_sin[1] * angle_sin[2];
        rot(2, 1) = angle_cos[0] * angle_sin[2]
                + angle_sin[0] * angle_sin[1] * angle_cos[2];
        rot(0, 2) = angle_sin[1];
        rot(1, 2) = -angle_cos[1] * angle_sin[2];
        rot(2, 2) = angle_cos[1] * angle_cos[2];

        // rotate the hydrogens and put the coordinates
        // we do not want to consider if the individual hydrogens
        // need to be replaced, because if h[0] has a correct bond
        // length, but is placed at the position where we predict h[1]
        // this leads to crashes later on. If one was not good, we 
        // replace both
        sys.mol(m).pos(h[0]) = sys.mol(m).pos(a) + rot*v1;
        count++;
        sys.mol(m).pos(h[1]) = sys.mol(m).pos(a) + rot*v2;
        count++;
      }
      break;
    }
    case(6):
    {
      // nh4+, simliar to case 5
      double bond1 = find_bond(sys, gff, m, Bond(a, h[0], false), 0.1);
      double bond2 = find_bond(sys, gff, m, Bond(a, h[1], false), 0.1);
      double bond3 = find_bond(sys, gff, m, Bond(a, h[2], false), 0.1);
      double bond4 = find_bond(sys, gff, m, Bond(a, h[3], false), 0.1);
      Vec v01 = sys.mol(m).pos(a) - sys.mol(m).pos(h[0]);
      Vec v02 = sys.mol(m).pos(a) - sys.mol(m).pos(h[1]);
      Vec v03 = sys.mol(m).pos(a) - sys.mol(m).pos(h[2]);
      Vec v04 = sys.mol(m).pos(a) - sys.mol(m).pos(h[3]);
      if (fabs(v01.abs() - bond1) / bond1 > eps ||
              fabs(v02.abs() - bond2) / bond2 > eps ||
              fabs(v03.abs() - bond3) / bond3 > eps ||
              fabs(v04.abs() - bond4) / bond4 > eps) {

        // no angle search, theta is 109.5
        // phi is 0, 120, 240
        double angle = M_PI / 180.0 * 109.5;
        Vec v1(0.0, 0.0, bond1);
        Vec v2(bond2 * sin(angle) * cos(0.0),
                bond2 * sin(angle) * sin(0.0),
                bond2 * cos(angle));
        Vec v3(bond3 * sin(angle) * cos(M_PI / 180.0 * 120.0),
                bond3 * sin(angle) * sin(M_PI / 180.0 * 120.0),
                bond3 * cos(angle));
        Vec v4(bond4 * sin(angle) * cos(M_PI / 180.0 * 240.0),
                bond4 * sin(angle) * sin(M_PI / 180.0 * 240.0),
                bond4 * cos(angle));

        //get three random numbers for the angles
        //calculate sin and cos of these angles

        Vec angle_cos, angle_sin;
        for (int i = 0; i < 3; i++) {
          double ang = 2.0 * M_PI / RAND_MAX * double(rand());
          angle_cos[i] = cos(ang);
          angle_sin[i] = sin(ang);
        }

        gmath::Matrix rot(3, 3);
        rot(0, 0) = angle_cos[0] * angle_cos[1];
        rot(1, 0) = angle_sin[0] * angle_cos[2]
                + angle_cos[0] * angle_sin[1] * angle_sin[2];
        rot(2, 0) = angle_sin[0] * angle_sin[2]
                - angle_cos[0] * angle_sin[1] * angle_cos[2];
        rot(0, 1) = -angle_sin[0] * angle_cos[1];
        rot(1, 1) = angle_cos[0] * angle_cos[2]
                - angle_sin[0] * angle_sin[1] * angle_sin[2];
        rot(2, 1) = angle_cos[0] * angle_sin[2]
                + angle_sin[0] * angle_sin[1] * angle_cos[2];
        rot(0, 2) = angle_sin[1];
        rot(1, 2) = -angle_cos[1] * angle_sin[2];
        rot(2, 2) = angle_cos[1] * angle_cos[2];



        // rotate the hydrogens and put the coordinates
        // we do not want to consider if the individual hydrogens
        // need to be replaced, because if h[0] has a correct bond
        // length, but is placed at the position where we predict h[1]
        // this leads to crashes later on. If one was not good, we 
        // replace all four
        sys.mol(m).pos(h[0]) = sys.mol(m).pos(a) + rot*v1;
        count++;
        sys.mol(m).pos(h[1]) = sys.mol(m).pos(a) + rot*v2;
        count++;
        sys.mol(m).pos(h[2]) = sys.mol(m).pos(a) + rot*v3;
        count++;
        sys.mol(m).pos(h[3]) = sys.mol(m).pos(a) + rot*v4;
        count++;
      }
      break;
    }
    case(7):
    {
      // charged NH, connected to 3 NH atoms.
      double bond = find_bond(sys, gff, m, Bond(a, h[0], false), 0.1);
      Vec v0 = sys.mol(m).pos(a) - sys.mol(m).pos(h[0]);
      if (fabs(v0.abs() - bond) / bond > eps) {

        Vec v1 = sys.mol(m).pos(nh[0]) - sys.mol(m).pos(a);
        Vec v2 = sys.mol(m).pos(nh[1]) - sys.mol(m).pos(a);
        Vec v3 = sys.mol(m).pos(nh[2]) - sys.mol(m).pos(a);
        Vec v4 = (v1 + v2 + v3).normalize();
        sys.mol(m).pos(h[0]) = sys.mol(m).pos(a) - bond*v4;
        count++;
      }
      break;
    }
    case(8):
    {
      // charged NH2, connected to 2 NH atoms
      double bond1 = find_bond(sys, gff, m, Bond(a, h[0], false), 0.1);
      double bond2 = find_bond(sys, gff, m, Bond(a, h[1], false), 0.1);
      Vec v0 = sys.mol(m).pos(a) - sys.mol(m).pos(h[0]);
      Vec v1 = sys.mol(m).pos(a) - sys.mol(m).pos(h[1]);
      if (fabs(v0.abs() - bond1) / bond1 > eps ||
              fabs(v1.abs() - bond2) / bond2 > eps) {

        Vec v2 = sys.mol(m).pos(nh[0]) - sys.mol(m).pos(a);
        Vec v3 = sys.mol(m).pos(nh[1]) - sys.mol(m).pos(a);
        Vec v4 = -(v2 + v3).normalize();
        Vec v5 = v2.cross(v3).normalize();
        double angle = M_PI / 180.0 * find_angle(sys, gff, m,
                Angle(h[0], a, h[1], false), 109.5);

        // we do not want to consider if the individual hydrogens
        // need to be replaced, because if h[0] has a correct bond
        // length, but is placed at the position where we predict h[1]
        // this leads to crashes later on. If one was not good, we 
        // replace both
        sys.mol(m).pos(h[0]) = sys.mol(m).pos(a) +
                bond1 * sin(0.5 * angle) * v5 +
                bond1 * cos(0.5 * angle) * v4;
        sys.mol(m).pos(h[1]) = sys.mol(m).pos(a) -
                bond2 * sin(0.5 * angle) * v5 +
                bond2 * cos(0.5 * angle) * v4;
        count++;
      }
      break;
    }


  }
  return count;
}

int utils::generate_hcoordinates(gcore::System &sys, gcore::GromosForceField &gff, int m, int a, double eps) {

        // divide into hydrogens and non-hydrogens
        vector<int> h;
        vector<int> nh;
        get_h_nh_neighbours(sys, gff, m, a, h, nh);

        // only continue if we have hydrogens
        int numH = h.size();
        int numNH = nh.size();
        int geom = get_geometry(numH, numNH);
        if (numH && !geom) {
          ostringstream os;
          os << "Unexpected geometry for hydrogen bound to atom: "
                  << m + 1 << ":" << a + 1 << endl;
          throw (gromos::Exception("gch", os.str()));
        }

        // we have to have a geometry (this means that there are hydrogens)
        // and a should not be a hydrogen itself. (in the case of H2O we have
        // e.g. H-H bonds. These are only treated via the oxygen.
        int r=0;
        if (geom) {
          r = generate_hcoordinates(sys, gff, m, a, h, nh, geom, eps);
        }
        return r;
}

double find_bond(System &sys, GromosForceField &gff, int m, Bond b, double guess) {
  BondIterator bi(sys.mol(m).topology());
  double value = 0.0;
  for (; bi; ++bi) {
    if (bi()[0] == b[0] && bi()[1] == b[1]) {
      value = gff.bondType(bi().type()).b0();
      break;
    }
  }
  if (value != 0.0)
    return value;
  else
    return guess;
}

double find_angle(System &sys, GromosForceField &gff, int m, Angle a, double guess) {
  AngleIterator ai(sys.mol(m).topology());
  double value = 0.0;

  for (; ai; ++ai) {
    if (ai()[0] == a[0] && ai()[1] == a[1] && ai()[2] == a[2]) {
      value = gff.angleType(ai().type()).t0();
      break;
    }
  }
  if (value != 0.0)
    return value;
  else
    return guess;
}

int find_dihedral(System &sys, int m, int i, int j, std::vector<int> h) {

  if (j < i) {
    int t = i;
    i = j;
    j = t;
  }

  DihedralIterator di(sys.mol(m).topology());
  int fourth = -1;

  for (; di; ++di) {
    if (di()[1] == i && di()[2] == j) {
      for (unsigned int k = 0; k < h.size(); k++) {
        if (di()[0] == h[k]) fourth = di()[3];
        if (di()[3] == h[k]) fourth = di()[0];
      }
      if (fourth != -1) break;
    }
  }
  if (fourth == -1) {

    // cannot find a dihedral, then just take the first atom bonded to i, which is not a hydrogen
    utils::Neighbours n(sys, m, i);
    for (unsigned int k = 0; k < n.size(); k++) {
      int hydrogen = 0;
      for (unsigned int l = 0; l < h.size(); l++)
        if (n[k] == h[l]) hydrogen = 1;
      if (!hydrogen && n[k] != j) fourth = n[k];
    }
    if (fourth == -1) {
      utils::Neighbours n(sys, m, j);
      for (unsigned int k = 0; k < n.size(); k++) {
        int hydrogen = 0;
        for (unsigned int l = 0; l < h.size(); l++)
          if (n[k] == h[l]) hydrogen = 1;
        if (!hydrogen && n[k] != i) fourth = n[k];
      }
    }
    if (fourth == -1)
      throw (gromos::Exception("find_dihedral", "undefined position"));
  }

  return fourth;

}

void utils::get_h_nh_neighbours(System &sys, GromosForceField &gff, int m, int a, std::vector<int> &h, std::vector<int> &nh){

    // find the neighbours of this atom
    utils::Neighbours n(sys, m, a);
    int dummyType = gff.dummyAtomType();
        for (unsigned int i = 0; i < n.size(); i++) {
          if (sys.mol(m).topology().atom(n[i]).isH()) {
            h.push_back(n[i]);
          } else {
            // only add it if is not a dummy atom
            if (dummyType != -1 &&
                    sys.mol(m).topology().atom(n[i]).iac() != dummyType)
              nh.push_back(n[i]);
          }
        }
}

int utils::get_geometry(int numH, int numNH) {
        int geom = 0;
        if (numH == 1 && numNH == 1) geom = 1;
        if (numH == 1 && numNH == 2) geom = 2;
        if (numH == 2 && numNH == 1) geom = 3;
        if (numH == 3 && numNH == 1) geom = 4;
        // crystallographic water
        if (numH == 2 && numNH == 0) geom = 5;
        // nh4+
        if (numH == 4 && numNH == 0) geom = 6;
        if (numH == 1 && numNH == 3) geom = 7;
        if (numH == 2 && numNH == 2) geom = 8;
        
        return geom;
}
