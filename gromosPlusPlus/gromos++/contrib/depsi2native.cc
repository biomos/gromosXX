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

/**
 * program depsi2native
 */

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/System.h"
#include "../src/gromos/Exception.h"
#include "../src/gmath/Vec.h"
#include "../src/args/BoundaryParser.h"
#include "../src/bound/Boundary.h"
#include "../src/gio/InG96.h"
#include "../src/utils/groTime.h"
#include "../src/gcore/Box.h"

const double NH_dist = 0.1;

using namespace std;
using namespace utils;
using namespace args;
using namespace gio;
using namespace gcore;
//using namespace gmath;
using namespace bound;

bool exchange(int n, vector<int> v);
Vec get_pos_CNC(AtomSpecifier & plane, Boundary *pbc, System sys);
Vec get_pos_OCO(AtomSpecifier & plane, Boundary *pbc, System sys);
Vec get_pos_M(AtomSpecifier & plane, Boundary *pbc, System sys);
void print_residue_coords(AtomSpecifier residue, Vec H, vector<int> resnum);

int main(int argc, char** argv) {
  
  try {

    Argument_List knowns;
    knowns << "topo" << "res" << "H_pos" << "pbc" << "time" << "traj";

    string usage = "# " + string(argv[0]);
    usage += "\n\t@topo      <molecular topology file of depsi-peptide>\n";
    usage += "\t@res       <residue number where an exchange OA -> NH will be performed>\n";
    usage += "\t@H_pos     <plane where H lies in: CNC, OCO or M (in between)>\n";
    usage += "\t@pbc       <boundary type>\n";
    usage += "\t[@time      <start time and time step>]\n";
    usage += "\t@traj      <trajectory file of depsi-peptide simulation>\n";

    Arguments args(argc, argv, knowns, usage);

    // read the topology and prepare the system
    args.check("topo", 1);
    InTopology it(args["topo"]);
    System sys = it.system();

    // get the periodic boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);

    // add all atoms to an AtomSpecifier atoms
    AtomSpecifier atoms(sys), residue(sys), prev_residue(sys), plane(sys);
    atoms.addSpecifier("a:a");

    // read the residue numbers where the OA -> NH change should be done
    vector<int> resnum;
    for (Arguments::const_iterator iter = args.lower_bound("res");
            iter != args.upper_bound("res"); ++iter) {
      resnum.push_back(atoi((iter->second).c_str()));
    }

    // loop over all trajectories
    InG96 ic;
    Time time(args);
    Vec H;
    for (Arguments::const_iterator iter = args.lower_bound("traj");
            iter != args.upper_bound("traj"); ++iter) {

      // open trajectory file
      ic.open((iter->second).c_str());
      ic.select("ALL");

      // write the header of the generated trajectory file
      cout << "TITLE" << endl;
      cout << "  GromosXX" << endl;
      cout << "  trajectory with added hydrogens (OA -> NH)" << endl;
      cout << "  the positions of the hydrogens were generated along method "
              << args["H_pos"] << endl;
      cout << "END" << endl;

      // loop over single trajectory
      while (!ic.eof()) {
        ic >> sys;// >> time;

        // print the time block
        cout << "TIMESTEP" << endl;
        cout << setw(15) << time.steps() << setw(15) << time.time() << endl;
        cout << "END" << endl;
        cout << "POSITIONRED" << endl;

        // loop over all atoms
        int resn = 0;
        for (int a = 0; a < atoms.size(); ++a) {
          // get the atoms of the actual residue
          if (resn == atoms.resnum(a)) {
            residue.addAtom(atoms.mol(a), atoms.atom(a));
          } else {
            if (exchange(resn, resnum)) {
              // all atoms of the actual residues are stored in AtomSpecifier residue
              // now get the atoms defining the plane(s):
              // C and O from the previous residue
              // OA and CA from the actual residue
              for (int r = 0; r < prev_residue.size(); ++r) {
                if (prev_residue.name(r) == "O" || prev_residue.name(r) == "C") {
                  plane.addAtom(prev_residue.mol(r), prev_residue.atom(r));
                }
              }
              for (int r = 0; r < residue.size(); ++r) {
                if (residue.name(r) == "OA" || residue.name(r) == "CA") {
                  plane.addAtom(residue.mol(r), residue.atom(r));
                }
              }
              if (plane.size() != 4) {
                ostringstream msg;
                msg << "Wrong number of O, C, OA and CA atoms in residue " << resn + 1
                        << ", the correct number in @H_pos OCO mode should be 4!";
                throw gromos::Exception(argv[0], msg.str());
              }
              // calculate the position of the added H atom
              if (args["H_pos"] == "CNC") {
                H = get_pos_CNC(plane, pbc, sys);
              } else if (args["H_pos"] == "OCO") {
                H = get_pos_OCO(plane, pbc, sys);
              } else if (args["H_pos"] == "M") {
                H = get_pos_M(plane, pbc, sys);
              } else {
                ostringstream msg;
                msg << "Method @H_pos " << args["H_pos"] << " unknown!";
                throw gromos::Exception(argv[0], msg.str());
              }
            }
            // print the frame
            print_residue_coords(residue, H, resnum);
            // prepare for the next residue
            plane.clear();
            prev_residue = residue;
            residue.clear();
            residue.addAtom(atoms.mol(a), atoms.atom(a));
            ++resn;
          }
        }
        // the last residue has not been printed yet
        print_residue_coords(residue, H, resnum);
        // prepare for the next frame
        plane.clear();
        prev_residue = residue;
        residue.clear();
        // and the box information shall be printed too
        cout << "END" << endl;
        cout << "GENBOX" << endl;
        cout.precision(9);
        cout << sys.box().boxformat() << endl;
        cout << fixed << setw(15) << sys.box().K().abs() << setw(15)
                << setw(15) << sys.box().L().abs() << setw(15) << sys.box().M().abs() << endl
                << setw(15) << sys.box().alpha() << setw(15) << sys.box().beta() << setw(15) << sys.box().gamma() << endl
                << setw(15) << 0.0 << setw(15) << 0.0 << setw(15) << 0.0 << endl
                << setw(15) << sys.box().X() << setw(15) << sys.box().Y() << setw(15) << sys.box().Z() << endl  
                << "END" << endl;
      }
      ic.close();
    }

  } catch (const gromos::Exception & e) {
    cerr << e.what() << endl;
    exit(1);
  }

  return 0;
} /* end of main function */

bool exchange(int n, std::vector<int> v) {
  for (unsigned int i = 0; i < v.size(); ++i) {
    if (n == v[i] - 1) return true;
  }
  return false;
}

Vec get_pos_CNC(AtomSpecifier & plane, Boundary *pbc, System sys) {
  Vec O, C, OA, CA;
  // the center position, the other positions will be nearest images to this position
  for (int i = 0; i < plane.size(); ++i) {
    if (plane.name(i) == "OA") OA = plane.pos(i);
  }
  // get the other positions (C, C) as nearest images to OA
  for (int i = 0; i < plane.size(); ++i) {
    if (plane.name(i) == "O") O = pbc->nearestImage(OA, plane.pos(i), sys.box());
    if (plane.name(i) == "C") C = pbc->nearestImage(OA, plane.pos(i), sys.box());
    if (plane.name(i) == "CA") CA = pbc->nearestImage(OA, plane.pos(i), sys.box());
  }
  // Calculate the position of the H atom
  Vec H = ((OA - CA).normalize() - (C - OA).normalize()).normalize();
  H = OA + H * NH_dist;

  return H;
}

Vec get_pos_M(AtomSpecifier & plane, Boundary *pbc, System sys) {

  // IDEA:
  // get the two normal vectors of the plane => line = 0 + s * v with v the
  // direction vector of the line, since we know the NH distance the problem
  // is solved then ...

  Vec O, C, OA, CA;
  // the center position, the other positions will be nearest images to this position
  for (int i = 0; i < plane.size(); ++i) {
    if (plane.name(i) == "OA") OA = plane.pos(i);
  }
  // get the other positions (C, C) as nearest images to OA
  for (int i = 0; i < plane.size(); ++i) {
    if (plane.name(i) == "O") O = pbc->nearestImage(OA, plane.pos(i), sys.box());
    if (plane.name(i) == "C") C = pbc->nearestImage(OA, plane.pos(i), sys.box());
    if (plane.name(i) == "CA") CA = pbc->nearestImage(OA, plane.pos(i), sys.box());
  }
  // Calculate the position of the H atom
  Vec H = (OA - CA).normalize() - (C - OA).normalize();
  Vec e1 = (C - OA).normalize();
  // constract a vector exactly in the middle of the O-C-OA and C-OA-CA plane
  Vec v1 = ((O - C) - ((O - C).dot(e1) * e1)).normalize();
  Vec v2 = ((CA - OA) - ((CA - OA).dot(e1) * e1)).normalize();
  Vec e2 = (v1 + v2).normalize();
  // the normal vector of the two planes the H must be in
  Vec n1 = e1.cross(e2);
  Vec n2 = (C - OA).cross(CA - OA);
  Vec n3 = H.cross(n2);
  // here the solution of the szstem of the two equation follows
  //   (I) n1 * (x, y, z) = 0
  //  (II) n3 * (x, y, z) = 0
  double a1 = n1[0], b1 = n1[1], c1 = n1[2];
  double a2 = n3[0], b2 = n3[1], c2 = n3[2];
  double x = 1;
  double z = (a1 * b2 - a2 * b1) / (b1 * c2 - b2 * c1);
  double y = -z * (a1 * c2 - a2 * c1) / (a1 * b2 - a2 * b1);
  Vec h(x, y, z);
  H = OA + (h.normalize() * NH_dist);
  Vec HH = OA - (h.normalize() * NH_dist);
  // return the solution in the right direction
  if ((HH - C).abs() > (H - C).abs()) {
    return HH;
  }
  return H;
}

Vec get_pos_OCO(AtomSpecifier & plane, Boundary *pbc, System sys) {

  // IDEA:
  // get the two normal vectors of the plane => line = 0 + s * v with v the
  // direction vector of the line, since we know the NH distance the problem
  // is solved then ...

  Vec O, C, OA, CA;
  // the center position, the other positions will be nearest images to this position
  for (int i = 0; i < plane.size(); ++i) {
    if (plane.name(i) == "OA") OA = plane.pos(i);
  }
  // get the other positions (C, C) as nearest images to OA
  for (int i = 0; i < plane.size(); ++i) {
    if (plane.name(i) == "O") O = pbc->nearestImage(OA, plane.pos(i), sys.box());
    if (plane.name(i) == "C") C = pbc->nearestImage(OA, plane.pos(i), sys.box());
    if (plane.name(i) == "CA") CA = pbc->nearestImage(OA, plane.pos(i), sys.box());
  }
  // Calculate the position of the H atom
  Vec n1 = (O - C).cross(OA - C);
  Vec H = (OA - CA).normalize() - (C - OA).normalize();
  Vec n2 = (C - OA).cross(CA - OA);
  Vec n3 = H.cross(n2);
  // here the solution of the szstem of the two equation follows
  //   (I) n1 * (x, y, z) = 0
  //  (II) n3 * (x, y, z) = 0
  double a1 = n1[0], b1 = n1[1], c1 = n1[2];
  double a2 = n3[0], b2 = n3[1], c2 = n3[2];
  double x = 1;
  double z = (a1 * b2 - a2 * b1) / (b1 * c2 - b2 * c1);
  double y = -z * (a1 * c2 - a2 * c1) / (a1 * b2 - a2 * b1);
  Vec h(x, y, z);
  H = OA + (h.normalize() * NH_dist);
  Vec HH = OA - (h.normalize() * NH_dist);
  // return the solution in the right direction
  if ((HH - C).abs() > (H - C).abs()) {
    return HH;
  }
  return H;
}

void print_residue_coords(AtomSpecifier residue, Vec H, vector<int> resnum) {
  static int count = 0;
  // set the counter to zero if a new frame begins
  if (residue.resnum(0) == 0) {
    count = 0;
    ;
  }
  cout.precision(9);
  for (int i = 0; i < residue.size(); ++i) {
    cout << fixed << setw(15) << residue.pos(i)[0] << setw(15) << residue.pos(i)[1]
            << setw(15) << residue.pos(i)[2] << endl;
    count++;
    if (count % 10 == 0) {
      cout << "#" << setw(10) << count << endl;
    }
    if (residue.name(i) == "OA" && exchange(residue.resnum(i), resnum)) {
      cout << fixed << setw(15) << H[0] << setw(15) << H[1] << setw(15) << H[2] << endl;
      count++;
      if (count % 10 == 0) {
        cout << "#" << setw(10) << count << endl;
      }
    }
  }
}
