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

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/Boundary.h"
#include "../src/gcore/System.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gromos/Exception.h"
#include "../src/utils/AtomicRadii.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"

/* Current status:
- testing on a set of known spheres shows good agreement
- testing on a configuration of lysozyme gave a 4x larger value than the old sasa program
- what could be the reason? inclusions?
- in a second phase: think what to do with inclusions (discuss at GROMOS meeting?)
*/

//const double pi = gmath::physConst.get_pi();
const double pi = 3.1415279;
// here all the stuff used for this program only is defined
namespace sasa {
  
  // index and type (begin=0 or end=1) for intersection point 
  class iP {
  public:
    int index;
    int type;
    iP(int index_, int type_);
    iP() {};
  };

  class vec2 {
  protected:
    double pos[2];
  public:
    vec2() {};
    vec2(double x, double y);
    void setPos(double x, double y);
    double get_x();
    double get_y();
    double norm2();
    vec2 operator+(vec2 v);
    vec2 operator-(vec2 v);
  };
  
  // vector to intersection point pair
  class pairVec {
  public:
    vec2 begin, end;
    double vb, ve;
    pairVec(vec2 v1, vec2 v2, double vb_, double ve_);
  };
  
  // a sphere to represent an atom with radius r
  class sphere {
  protected:
    double r;
    Vec pos;
  public:
    sphere(double r_, Vec pos_);
    Vec get_pos();
    double get_radius();
  };

  class arc {
  protected:
    vec2 pos;
    double r; // radius of the arc
    double R; // radius of the vdW shell of the atom
    double dist; // distance between center of sphere and z-plane
    map<double, iP> mapPoints;
    vector<pairVec> overlapPoints;
    double len;
  public:
    arc(double x_, double y_, double r_, double len_, double R_, double dist_);
    double get_r();
    double get_R();
    double get_dist();
    double get_x();
    double get_y();
    vec2 get_pos();
    double get_len();
    void set_len(double l);
    void overlap(vector<arc>::iterator a1, vector<arc>::iterator a2);
    double calculate_angle(vec2 v, double r);
    void remove_intersections();
  };

}

using namespace args;
using namespace bound;
using namespace gcore;
using namespace gio;
using namespace sasa;
using namespace std;
using namespace utils;


int main(int argc, char ** argv) {

  try {

    if(true){
      std::cout << "\nThis program is under development\n"
                << "please do not use me\n\n";
      exit(0);
    }
    // the (known) arguments
    Argument_List knowns;
    knowns << "topo" << "atoms" << "time" << "probe" << "zslice" << "pbc" << "traj";

    string usage = "# " + string(argv[0]);
    usage += "\n\t@topo      <molecular topology file>\n";
    usage += "\t@atoms    <solut atoms to be considered to calculate the SASA>\n";
    usage += "\t[@probe    <the IAC and LJ radius of the solvent/probe; first solvent atom is taken if not specified>]\n";
    usage += "\t[@time     <time and dt>]\n";
    usage += "\t@pbc      <periodic boundary and gathering>\n";
    usage += "\t[@zslice   <distance between the Z-slices (default: 0.005)>]\n";
    usage += "\t@traj     <trajectory files>\n";

    Arguments args(argc, argv, knowns, usage);

    // read in the topology to build the system
    args.check("topo", 1);
    InTopology it(args["topo"]);
    System sys(it.system());

    // get the atom/LJ radius of the solvent
    std::vector<double> solvent = args.getValues<double>("probe", 2, false, 
            Arguments::Default<double>() << sys.sol(0).topology().atom(0).iac() 
                                         << sys.sol(0).topology().atom(0).radius());
    double solviac = solvent[0];
    double solvrad = solvent[1];

    // calculate the van der Waals radii
    compute_atomic_radii_vdw(solviac, solvrad, sys, it.forceField());

    // get the list of atoms to be considered for the SASA calculation
    AtomSpecifier atoms(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("atoms");
      Arguments::const_iterator to = args.upper_bound("atoms");
      for (; iter != to; iter++) {
        atoms.addSpecifier(iter->second.c_str());
      }
    }
    if (atoms.size() < 1) {
      stringstream msg;
      throw gromos::Exception("sasa", "no atoms found for given atom specifier (@atoms)");
    }

    // is there a trajectory file?
    if (args.count("traj") < 1) {
      throw gromos::Exception("sasa", "no trajectory file specified (@traj)");
    }

    // get the boundary and gather method
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, pbc->refSys(), args);

    // get the simulation time
    Time time(args);

    // get the distance between the x/y-planes
    double dz = args.getValue<double>("zslice", false, 0.005);
    
    // the output
    cout << "# Time   SASA" << endl;

    // loop over the trajectory file
    Arguments::const_iterator iter = args.lower_bound("traj");
    Arguments::const_iterator to = args.upper_bound("traj");
    for (; iter != to; ++iter) {
      InG96 ic;
      ic.open((iter->second).c_str());
      // loop over the different configurations/frames of the trajectory file
      while (!ic.eof()) {

        // read the current frame/configuration
        ic >> sys >> time;

        // gather the current configurations
        (*pbc.*gathmethod)();

        // build a sphere for each (solute) atom to be considered
        vector<sphere> spheres;
        for (int a = 0; a < atoms.size(); ++a) {
          // take only heavy atoms!
          if (!sys.mol(atoms.mol(a)).topology().atom(atoms.atom(a)).isH()) {
            double r = atoms.radius(a);
            Vec pos = atoms.pos(a);
            spheres.push_back(sphere(r, pos));
          }
        }

        // get z_min and z_max (to span the x/y-planes)
        double z_min = spheres[0].get_pos()[2] -
                spheres[0].get_radius() - solvrad;
        double z_max = spheres[0].get_pos()[2] +
                spheres[0].get_radius() + solvrad;
        for (unsigned int s = 0; s < spheres.size(); ++s) {
          double z = spheres[s].get_pos()[2];
          if ((z - spheres[s].get_radius() - solvrad) < z_min) {
            z_min = z - spheres[s].get_radius() - solvrad;
          } else if ((z + spheres[s].get_radius() + solvrad) > z_max) {
            z_max = z + spheres[s].get_radius() + solvrad;
          }
        }

        // overlap the spheres with the N x/y-planes
        // first get the number of planes needed, N
        int N = floor((z_max - z_min) / dz);
        if (N == 0) {
          N = 1;
        }
        
        // get a vector containing the circles/arcs per plane
        vector<vector<arc> > arcs(N);
        // this is the solvent accessible surface area;
        double area = 0.0;
        // loop over planes
        for (int n = 0; n < N; ++n) {
          // the z-position of the current plane
          double z0 = z_min + (n + 1) * dz; // we don't start at z_min since there
                                            // is by definition no sphere-plane overlap
          double length = 0.0;
          // loop over all spheres and add a circle/arc if the sphere cuts the current plane
          for (unsigned int s = 0; s < spheres.size(); ++s) {
            double R = spheres[s].get_radius() + solvrad; // sphere radius (including solvent/probe radius)
            double z_sphere = spheres[s].get_pos()[2]; // z-position of centre of sphere
            double dist = abs(z0 - z_sphere); // distance of centre of sphere to the current plane
            if (dist < R) { // the sphere cuts the current plane
              double x = spheres[s].get_pos()[0]; // x-position of circle/arc
              double y = spheres[s].get_pos()[1]; // y-position of circle/arc
              double r = R * sin(acos(dist / R)); // radius of circle/arc
              double len = 2 * pi * r; // the length of the arc, actually the whole circle
              arc a(x, y, r, len, R, dist);
              arcs[n].push_back(a);
            }
          } // end of loop over spheres

          // now overlap the circles within one plane and remove intersections
          // loop over circles/arcs within a plane
          vector<arc>::iterator it1 = arcs[n].begin();
          for (int a1 = 0; a1 < ((int) arcs[n].size() - 1); ++a1, ++it1) {
            vector<arc>::iterator it2 = it1;
            ++it2;
            //cerr << "circle " << a1 << ": " << it1->get_x() << " " << it1->get_y() << endl;
            for (unsigned int a2 = a1 + 1; a2 < arcs[n].size(); ++a2, ++it2) {
              it1->overlap(it1, it2);
            } // end of loop over (second) circles
            // now remove intersections and add up the length
            if (it1->get_len() > 0.0) 
              it1->remove_intersections();
            double factor = it1->get_R() / sqrt(it1->get_R()*it1->get_R() 
                            - it1->get_dist()*it1->get_dist());
            double zres = dz/2.0;
            if (it1->get_R() - it1->get_dist() < zres)
              zres += it1->get_R() - it1->get_dist();
            else 
              zres += zres;
            length += it1->get_len() * factor * zres;
          } // end of loop over (first) circles
          // add the last circle - if there is one
          if (it1 != arcs[n].end()) {
            double zres = dz / 2.0;
            double factor = it1->get_R() / sqrt(it1->get_R() * it1->get_R()
                    - it1->get_dist() * it1->get_dist());
            if (it1->get_R() - it1->get_dist() < zres)
              zres += it1->get_R() - it1->get_dist();
            else
              zres += zres;
            length += it1->get_len() * factor * zres;
          }
          area += length;
        } // end of loop over planes
        
        // print out the sasa
        cout << time.time() << "     " << area << endl;

      } // end of loop over configurations/frames
    } // end of loop over trajectory files
    

  } catch (const gromos::Exception &e) {
    // quit with an error message
    cerr << e.what() << endl;
    exit(1);
  }

  return 0;
}

// -----------------------------------------------------------------------------
// ======================= Function, Class, ... Definitions ====================
//
// first all the stuff in the sasa namespace
namespace sasa {
  
  iP::iP(int index_, int type_) {
    index = index_;
    type = type_;
  }

  vec2::vec2(double x, double y) {
    pos[0] = x;
    pos[1] = y;
  }
  
  void vec2::setPos(double x, double y) {
    pos[0] = x;
    pos[1] = y;
  }
  
  double vec2::get_x() {
    return pos[0];
  }
  
  double vec2::get_y() {
    return pos[1];
  }
  
  double vec2::norm2() {
    return pos[0] * pos[0] + pos[1] * pos[1];
  }
  
  vec2 vec2::operator+(vec2 v) {
    vec2 r(v.get_x() + get_x(), v.get_y() + get_y());
    return r;
  }
  
  vec2 vec2::operator-(vec2 v) {
    v.setPos(v.get_x() - pos[0], v.get_y() - pos[1]);
    return v;
  }
  
  pairVec::pairVec(vec2 v1, vec2 v2, double vb_, double ve_) {
    begin = v1;
    end = v2;
    vb = vb_;
    ve = ve_;
  }
  
  sphere::sphere(double r_, Vec pos_) {
    r = r_;
    pos = pos_;
  }

  Vec sphere::get_pos() {
    return pos;
  }

  double sphere::get_radius() {
    return r;
  }

  arc::arc(double x_, double y_, double r_, double len_, double R_, double dist_) {
    pos.setPos(x_, y_);
    r = r_;
    len = len_;
    R = R_;
    dist = dist_;
  }

  double arc::get_r() {
    return r;
  }
  
  double arc::get_R() {
    return R;
  }
  
  double arc::get_dist() {
    return dist;
  }
  
  double arc::get_x() {
    return pos.get_x();
  }
  
  double arc::get_y() {
    return pos.get_y();
  }

  vec2 arc::get_pos() {
    return pos;
  }
  
  double arc::get_len() {
    return len;
  }
  
  void arc::set_len(double l) {
    len = l;
  }

  void arc::overlap(vector<arc>::iterator a1, vector<arc>::iterator a2) {
    double x1 = a1->get_x();
    double y1 = a1->get_y();
    double r1 = a1->get_r();
    double x2 = a2->get_x();
    double y2 = a2->get_y();
    double r2 = a2->get_r();
    // distance of the two centres
    double dx = x1 - x2;
    double dy = y1 - y2;
    double d = sqrt(dx * dx + dy * dy);
    // now treat the different cases of overlap and no overlap
    if (d >= r1 + r2) { // no overlap, circles too far away
      return;
    } else if (abs(r1 - r2) >= d) { // there is no overlap, one circle in the other
      if (r1 < r2) {
        a1->set_len(0.0);
      } else {
        a2->set_len(0.0);
      }
    } else { // there is an overlap
        // overlap the two circles:
      // r1^2 = (x - x1)^2 + (y - y1)^2     (I)
      // r2^2 = (x - x2)^2 + (y - y2)^2    (II)
      // (II) - (I) => y = mx + n with    (III)
      double inter_x1, inter_x2, inter_y1, inter_y2;
      if (x1 == x2) {
        double m = -dx / dy;
        double R = (r1*r1 - r2*r2) - (y1*y1 - y2*y2) - (x1*x1 - x2*x2);
        double n = R / (2 * -dy);
        // insertion of y = mx + n in (I) leads to x_{1,2}
        // calculation of x_{1,2} using the p-q-formula:
        double m2i = 1.0 / (1 + m * m);
        double p = (2 * m * (n - y1) - 2 * x1) * m2i;
        double q = ((n - y1)*(n - y1) - r1 * r1 + x1 * x1) * m2i;
        double factor2 = sqrt((p * p / 4.0) - q);
        double factor1 = -p / 2.0;
        inter_x1 = factor1 + factor2;
        inter_x2 = factor1 - factor2;
        // insert x_{1,2} into (III)
        inter_y1 = m * inter_x1 + n;
        inter_y2 = m * inter_x2 + n;
      } else {
        double m = -dy / dx;
        double R = (r1*r1 - r2*r2) - (x1*x1 - x2*x2) - (y1*y1 - y2*y2);
        double n = R / (2 * -dx);
        // insertion of y = mx + n in (I) leads to x_{1,2}
        // calculation of x_{1,2} using the p-q-formula:
        double m2i = 1.0 / (1 + m * m);
        double p = (2*m * (n - x1) - 2*y1) * m2i;
        double q = ((n - x1)*(n - x1) - r1*r1 + y1*y1) * m2i;
        double factor2 = sqrt((p * p / 4.0) - q);
        double factor1 = -p / 2.0;
        inter_y1 = factor1 + factor2;
        inter_y2 = factor1 - factor2;
        // insert x_{1,2} into (III)
        inter_x1 = m * inter_y1 + n;
        inter_x2 = m * inter_y2 + n;
      }
      
      // STORE INTERSECTION POINTS
      vec2 vb1(inter_x1 - x1, inter_y1 - y1);
      vec2 ve1(inter_x2 - x1, inter_y2 - y1);
      vec2 vb2(inter_x1 - x2, inter_y1 - y2);
      vec2 ve2(inter_x2 - x2, inter_y2 - y2);
      
      if (r1 > r2) { // circle 1 larger than circle 2
        if ((vb1.get_y() >= 0 && ve1.get_y() < 0) && (vb1.get_x() >= 0 && ve1.get_x() >= 0)) {
          // exchange intersection points
          vec2 tmp = vb1;
          vb1 = ve1;
          ve1 = tmp;
        }
      } else if (r2 > r1) { // circle 2 larger than circle 1
        if ((vb2.get_y() >= 0 && ve2.get_y() < 0) && (vb2.get_x() >= 0 && ve2.get_x() >= 0)) {
          // exchange intersection points
          vec2 tmp = vb2;
          vb2 = ve2;
          ve2 = tmp;
        }
      } else { // the circles have the same size
        if ((vb1.get_y() >= 0 && ve1.get_y() < 0) && (vb1.get_x() >= 0 && ve1.get_x() >= 0)) {
          // exchange intersection points
          vec2 tmp = vb1;
          vb1 = ve1;
          ve1 = tmp;
        }
        if ((vb2.get_y() >= 0 && ve2.get_y() < 0) && (vb2.get_x() >= 0 && ve2.get_x() >= 0)) {
          // exchange intersection points
          vec2 tmp = vb2;
          vb2 = ve2;
          ve2 = tmp;
        }
      }
      // first arc
      // calculate angles
      double angle_begin = calculate_angle(vb1, r1);
      double angle_end = calculate_angle(ve1, r1);
      int index = a1->overlapPoints.size();
      iP ip_b1(index, 0);
      iP ip_e1(index, 1);
      pairVec pv1(vb1, ve1, angle_begin, angle_end);
      a1->overlapPoints.push_back(pv1);
      a1->mapPoints[angle_begin] = ip_b1;
      a1->mapPoints[angle_end] = ip_e1;
      //cerr << "arc1: begin = " << angle_begin << ", end = " << angle_end << endl;
      
      // second arc
      // calculate angles
      angle_begin = calculate_angle(vb2, r2);
      angle_end = calculate_angle(ve2, r2);
      index = a2->overlapPoints.size();
      iP ip_b2(index, 0);
      iP ip_e2(index, 1);
      pairVec pv2(vb2, ve2, angle_begin, angle_end);
      a2->overlapPoints.push_back(pv2);
      a2->mapPoints[angle_begin] = ip_b2;
      a2->mapPoints[angle_end] = ip_e2;
      //cerr << "arc2: begin = " << angle_begin << ", end = " << angle_end << endl;
    }
  }
  
  double arc::calculate_angle(vec2 v, double r) {
    if (v.get_y() < 0) {
      return  2 * pi - acos(v.get_x() / r);
    } else {
      return acos(v.get_x() / r);
    }
  }

  void arc::remove_intersections() {
    map<double, iP>::const_iterator it = mapPoints.begin();
    map<double, iP>::const_iterator to = mapPoints.end();
    double start_angle, end_angle;
    double limit_angle = 2 * pi;
    while (it != to && it->first < limit_angle) {
      // check whether the first point is a begin or an end
      if (it->second.type == 1) { // it's an end --> angle_begin > angle_end
        start_angle = overlapPoints[it->second.index].vb;
        end_angle = it->first;
        limit_angle = start_angle;
      } else {
        start_angle = it->first;
        end_angle = overlapPoints[it->second.index].ve;
      }
      //cerr << "remove: start = " << start_angle << ", end = " << end_angle << ", limit = " << limit_angle << endl;
      // now search for the next point
      it++;
      while (it != to && it->first < end_angle && it->first <= limit_angle) {
        if (it->second.type == 0) { // it's a beginning
          if (overlapPoints[it->second.index].ve > end_angle) { // the end is outside the range
            end_angle = overlapPoints[it->second.index].ve;
          }
        } else { // it's an end
          if (overlapPoints[it->second.index].vb > end_angle) { // --> angle_begin > angle_end
            start_angle = overlapPoints[it->second.index].vb;
            limit_angle = start_angle;
          }
        }
        //cerr << "start = " << start_angle << ", end = " << end_angle << ", limit = " << limit_angle << endl;
        it++;
      }
      // we found a pair, remove the length between the vectors
      // calculate the angle
      vec2 vb = overlapPoints[mapPoints[start_angle].index].begin;
      vec2 ve = overlapPoints[mapPoints[end_angle].index].end;
      double final_angle = acos((vb.get_x()*ve.get_x() + vb.get_y()*ve.get_y()) / (r*r));
      if (vb.get_x() >=0 && ve.get_y() >= 0 && vb.get_y() >=0 && ve.get_y() < 0)
        final_angle = 2*pi - final_angle;
      len -= r*final_angle;
      //cerr << "final = " << final_angle << ", len = " << len << ", r*final = " << r*final_angle << endl;
    }
  }
}
