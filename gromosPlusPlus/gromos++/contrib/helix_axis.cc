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
 * @file helix_axis.cc
 * find axis and dimensions of an alpha helix
 * based on Kahn, Computers Chem, 13:185-189, 1989
 * note that only the simplest fitting method, from one end of the helix to the other
 * is implemented, and that currently you must pass only the coordinates of the
 * helix residues
 */
/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor helix_axis
 * @section helix_axis
 * @author @ref ja
 * @date 19.04.2010
 *
 * PROGRAM DESCRIPTION
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> [\@time</td><td>&lt;time and dt&gt;] </td></tr>
 * <tr><td> [\@timespec</td><td>&lt;timepoints to consider for the fit: ALL (default), EVERY or SPEC (if time-series)&gt;] </td></tr>
 * <tr><td> [\@timepts</td><td>&lt;timepoints to consider for the fit (if time-series and timespec EVERY or SPEC)&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
   helix_axis
     @topo       ex.top
     @pbc        r
     @time       0 0.02
     @traj       ex.trj
   @endverbatim

 * <hr>
 */

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/args/ReferenceParser.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Matrix.h"
#include "../src/gmath/Physics.h"
#include "../src/utils/groTime.h"
#include "../src/gromos/Exception.h"

using namespace args;
using namespace bound;
using namespace fit;
using namespace gcore;
using namespace gio;
using namespace gmath;
using namespace std;
using namespace utils;

// function for skipping time-points
bool use_this_frame(int i, string const & timespec, vector<int> const & timepts,
        unsigned int & timesComputed, bool & done);

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "timespec" << "timepts" << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo      <molecular topology file>\n";
  usage += "\t@pbc         <boundary type> [<gathermethod>]\n";
  usage += "\t[@time       <t> <dt>] (optional; will only print time series if given)\n";
  usage += "\t[@timespec   <timepoints at which to compute the SASA: ALL (default), EVERY or SPEC>]\n";
  usage += "\t[@timepts    <timepoints at which to compute the SASA>] (if timespec EVERY or SPEC)\n";
  usage += "\t@traj        <trajectory files>\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // first print a warning:
    cout << "# WARNING: this program is experimental. Please consult Jane before using\n" << endl;

    // define physical constants
    const double radian2degree = gmath::physConst.get_radian2degree();

    // read topology and initialise system
    InTopology it(args["topo"]);
    System sys(it.system());

    System refSys(it.system());

    // now parse boundary conditions for sys
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);
    
    // get time
    Time time(args);

    // parse timespec
    string timespec = "ALL";
    vector<int> timepts;
    if (args.count("timespec") > 0) {
      timespec = args["timespec"];
      if (timespec != "ALL" && timespec != "EVERY" && timespec != "SPEC")
        throw gromos::Exception("helix_axis",
              "timespec format " + timespec + " unknown.\n");
      if (timespec == "EVERY" || timespec == "SPEC") {
        for (Arguments::const_iterator it = args.lower_bound("timepts");
                it != args.upper_bound("timepts"); ++it) {
          int bla = atoi(it->second.c_str());
          timepts.push_back(bla);
        }
        if (timepts.size() == 0) {
          throw gromos::Exception("helix_axis",
                  "if you give EVERY or SPEC you have to use "
                  "@timepts as well");
        }
        if (timepts.size() != 1 && timespec == "EVERY") {
          throw gromos::Exception("helix_axis",
                  "if you give EVERY you have to give exactly"
                  " one number with @timepts");
        }
      }
    }

    // check we have enough residues (nres is true number)
    // note only one molecule allowed
    int nres = sys.mol(0).topology().numRes();
    if (nres < 4) {
      throw gromos::Exception("helix_axis",
              "Unable to define helix axis with fewer than four residues");
    } else if (nres < 6) {
      cout << "Warning: helix axis definition is shaky with less than 6 residues" << endl;
    }
    // define CA atom numbers
    int CA1, CA2, CA3, CAp, CApm1, CApp1;
    // find second and penultimate CA atom numbers (name == CA)
    int natoms = sys.mol(0).numAtoms();
    for (int ii = 0; ii < natoms; ++ii) {
      int resnum = sys.mol(0).topology().resNum(ii);
      string atomtype = sys.mol(0).topology().atom(ii).name();
      if (atomtype == "CA") {
        if (resnum == 0) {
          CA1 = ii;
        } else if (resnum == 1) {
          CA2 = ii;
        } else if (resnum == 2) {
          CA3 = ii;
        } else if (resnum == nres - 3) {
          CApm1 = ii;
        } else if (resnum == nres - 2) {
          CAp = ii;
        } else if (resnum == nres - 1) { // actually the last residue
          CApp1 = ii;
        }
      }
    }
    // check that CA atoms were found
    cout << "CA atoms " << CA1 << " " << CA2 << " " << CA3 << " " << CApm1 <<
            " " << CAp << " " << CApp1 << endl;

    // define origin
    Vec origin(0.0, 0.0, 0.0);
    // define z-axis
    Vec zaxis(0.0, 0.0, 1.0);

    // start at -1 to get times right
    int num_frames = -1;
    // number of time-points for which helix axis and dimensions have been calculated
    unsigned int times_computed = 0;
    // for SPEC: so that we stop trying when all requested timepoints are written
    bool done = false;

    // initialise output
    cout << "#\tTime\tRadius\t\tLength" << endl;

    // initialise averages
    double nvalues = 0.0; // duplicates above counters but needs to be double
    double rave = 0.0;
    double lave = 0.0;

    // loop over all trajectories
    if (args.count("traj") > 0) {
      for (Arguments::const_iterator
        iter = args.lower_bound("traj"),
              to = args.upper_bound("traj");
              iter != to; ++iter) {

        // open file
        InG96 ic;
        ic.open((iter->second).c_str());
        ic.select("SOLUTE");

        // loop over single trajectory
        while (!ic.eof()) {

          ic >> sys >> time;
          // check the input
          if (sys.hasBox == false && pbc->type() != 'v')
            throw gromos::Exception("helix_axis",
                  "If pbc != v the trajectory file should have a box block ");
          if (!sys.hasPos)
            throw gromos::Exception("helix_axis",
                  "Unable to read POSITION(RED) block from "
                  "trajectory file.");

          // gather
          (*pbc.*gathmethod)();

          // check whether to skip this frame or not
          num_frames++;
          if (use_this_frame(num_frames, timespec, timepts, times_computed, done)) {

            // vector from CA2 to origin
            Vec P1 = origin - sys.mol(0).pos(CA2);
            // vector from CA2 to CA1
            Vec A1 = sys.mol(0).pos(CA2) - sys.mol(0).pos(CA1);
            // vector from CA2 to CA3
            Vec B1 = sys.mol(0).pos(CA2) - sys.mol(0).pos(CA3);
            // bisector of angle A1B1
            Vec V1 = (A1 + B1).normalize();

            // vector from penultimate CA (CAp) to origin
            Vec P2 = origin - sys.mol(0).pos(CAp);
            // vector from CAp to CApm1 (penultimate minus 1)
            Vec A2 = sys.mol(0).pos(CAp) - sys.mol(0).pos(CApm1);
            // vector from CAp to CApp1 (penultimate plus 1)
            Vec B2 = sys.mol(0).pos(CAp) - sys.mol(0).pos(CApp1);
            // bisector of angle A2B2
            Vec V2 = (A2 + B2).normalize();

            // helix axis is cross-product of V1 and V2
            Vec H = V1.cross(V2).normalize();
            cout << "H " << H[0] << " " << H[1] << " " << H[2] << endl;

            // d is distance along axis between CA2 and CAp
            double d = (P2 - P1).dot(H);
            // check components again
            double dhabs2 = (d*H).abs2();
            double P2P1abs2 = (P2-P1).abs2();
            double P2P1V2abs = fabs((P2-P1).dot(V2));
            cout << "Radius components " << d << " " << dhabs2 << " " << P2P1abs2 << " " <<
                    P2P1V2abs << endl;
            // compute radius
            //double r = fabs(((d * H).abs2() - (P2 - P1).abs2()) /
            //(2. * (fabs((P1 - P2).dot(V2)))));
            // try with P2 - P1 in denominator
            double r = fabs(((d * H).abs2() - (P2 - P1).abs2()) /
                    (2. * (fabs((P2 - P1).dot(V2)))));
            rave += r;

            // find angle between H and z-axis (both have length 1.0)
            double angle = acos(H.dot(zaxis)) * radian2degree;
            // find vector perpendicular to both H and z-axis
            Vec perp = H.cross(zaxis);

            // find rotation matrix
            Matrix rotmat = PositionUtils::rotateAround(perp, angle);
            // rotate helix so helix-axis(H) is parallel to z-axis of coordinate system
            PositionUtils::rotate(&sys, rotmat);

            // compute length of helix as diff in z coordinates
            double maxz = PositionUtils::getmaxcoordinates(&sys, false)[2];
            double minz = PositionUtils::getmincoordinates(&sys, false)[2];
            double length = maxz - minz;
            lave += length;

            // write out time, radius (r) and length of helix (dir now == z)
            cout << time << "\t" << r << "\t" << length << endl;

            // increment nvalues
            nvalues += 1.0;

          }// use this frame
        }// frames loop
      }// trajs loop

      // write out averages
      cout << "#Averages:\t" << rave / nvalues << "\t" << lave / nvalues << endl;

    } else {
      throw gromos::Exception("helix_axis",
              "You need to give a trajectory file");
    }



  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

// function to decide whether to use this frame or not
bool use_this_frame(int i, std::string const & timespec, vector<int> const & timepts,
        unsigned int & timesComputed, bool & done) {
  if (timespec == "ALL") {
    ++timesComputed;
    return true;
  } else if (timespec == "EVERY" && i % timepts[0] == 0) {
    ++timesComputed;
    return true;
  } else if (timespec == "SPEC") {
    for (unsigned int j = 0; j < timepts.size(); ++j) {
      if (timepts[j] == i) {
        ++timesComputed;
        if (timesComputed == timepts.size())
          done = true;
        return true;
      } // compute
    } // times
  }
  return false;
}
