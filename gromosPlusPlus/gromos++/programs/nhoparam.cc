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
 * @file nhoparam.cc
 * calculates order parameters for N-H bonds
 */
/**
 * @page programs Program Documentation
 *
 * @anchor nhoparam
 * @section nhoparam calculates order parameters for N-H bonds in protein backbones
 * @author @ref mk @ref ns
 * @date 23.10.2008
 *
 * Program nhoparam calculates order parameters for a given set of nitrogen atoms
 * (@ref AtomSpecifier \@atoms).
 *
 * In a first step, the program determines the N-H bonds (of which @f$\mu@f$ 
 * is the unit vector) by the atomic masses of nitrogen and hydrogen. For 
 * secondary and tertiary amides the different N-H bonds are averaged. Then,
 *
 * @f[ S^2 = \frac{1}{2} \left[ 3 \sum_{i=1}^3 \sum_{j=1}^3  \left< \mu_{i}(t)\mu_{j}(t) \right>_{t}^2 - 1\right]@f]
 *
 * is applied in order to calculate the order parameter of the N-H bond after performing 
 * a least-square rotational fit. Fitting can be controlled using the \@ref and
 * \@atomsfit arguments. If \@ref is absent, the first frame of the trajectory is
 * taken as reference. \@atomsfit are the @ref AtomSpecifier "atoms" used for
 * fitting. If omitted, the nitrogen atoms are used.
 *
 * The running averaged and window averaged (using a window size of \@winframe)
 * order parameters are written to two seperate time series files (OPts.out, 
 * OPwints.out). Final results and statistics are written to the standard output.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gathermethod&gt;] </td></tr>
 * <tr><td>[\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;]</td></tr>
 * <tr><td> \@winframe</td><td>&lt;averaging window (number of frames)&gt; </td></tr>
 * <tr><td> \@atoms</td><td>&lt;nitrogen @ref AtomSpecifier "atom(s)" to calculate order parameter from&gt; </td></tr>
 * <tr><td> [\@atomsfit</td><td>&lt;@ref AtomSpecifier "atoms" to consider for fit&gt;] </td></tr>
 * <tr><td> [\@ref</td><td>&lt;reference coordinates (if absent, the first frame of \@traj is reference)&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  nhoparam
    @topo ex.top
    @pbc r
    @atoms 1:N
    @ref exref.coo
    @atomsfit 1:CA,C,N
    @time 0 0.1
    @winframe 10
    @traj ex.tr
    @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <string>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <set>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/args/ReferenceParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/utils/Neighbours.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Matrix.h"
#include "../src/gmath/Stat.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;
using namespace fit;

/**
 * calculates the S2 value from the sum of the tensor elements and
 * number of frames
 */
double S2_calc(const gmath::Matrix & sum, int num_frames);

// masses for checking of the N-H bonds.
static const double nitrogen_mass = 14.00670;
static const double hydrogen_mass = 1.00800;

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "traj" << "atoms" << "ref" << "pbc" << "atomsfit"
          << "time" << "winframe";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@atoms <nitrogen atom(s) to calculate order parameter from>\n";
  usage += "\t[@ref <reference coordinates>]\n";
  usage += "\t[@atomsfit <atoms to consider for fit>]\n";
  usage += "\t@time <time and dt>\n";
  usage += "\t@winframe <averaging window [# of frames]>\n";
  usage += "\t@traj <trajectory files>\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System refSys(it.system());

    System sys(refSys);

    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(refSys, args);
    // GatherParser
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

    // read reference coordinates...
    InG96 ic;
    if (args.count("ref") > 0)
      ic.open(args["ref"]);
    else
      if (args.count("traj") > 0)
      ic.open(args.lower_bound("traj")->second);

    ic >> refSys;
    ic.close();

    // this always goes wrong. check that there is a box block in the refSys
    if (refSys.hasBox == false && pbc->type() != 'v')
      throw gromos::Exception(argv[0],
            "If pbc != v you have to give a box block "
            "in the reference system as well.");
    // and that we managed to read coordinates from it
    if (!refSys.hasPos)
      throw gromos::Exception(argv[0],
            "Unable to read POSITION(RED) block from "
            "reference positions file.");

    // gather reference system
    (*pbc.*gathmethod)();

    delete pbc;

    // System for calculation
    //System sys(refSys);
    AtomSpecifier fitatoms(refSys);
    AtomSpecifier atoms(sys);


    //get nitrogen atoms
    {
      Arguments::const_iterator iter = args.lower_bound("atoms");
      Arguments::const_iterator to = args.upper_bound("atoms");

      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        atoms.addSpecifier(spec);
      }
    }
    if (atoms.size() == 0)
      throw gromos::Exception(argv[0], "No nitrogen atoms specified!");

    //try for fit atoms
    if (args.count("atomsfit") > 0) {
      Arguments::const_iterator iter = args.lower_bound("atomsfit");
      Arguments::const_iterator to = args.upper_bound("atomsfit");

      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        fitatoms.addSpecifier(spec);
      }
    } else {
      cout << "# @atoms atoms are taken as fit atoms." << endl;
      fitatoms = atoms;
    }
    RotationalFit * rf = NULL;
    if (fitatoms.size()) {
      Reference * ref = new Reference(&refSys);
      ref->addAtomSpecifier(fitatoms);
      rf = new RotationalFit(ref);
    }


    // Parse boundary conditions for sys
    pbc = BoundaryParser::boundary(sys, args);
    // create the rotational fit

    //see whether no atoms are supplied
    vector<set<int> > hydrogens(atoms.size(), set<int>());
    for (unsigned int i = 0; i < atoms.size(); ++i) {
      if (atoms.mass(i) != nitrogen_mass) {
        ostringstream msg;
        msg << "Atom " << atoms.toString(i) << " is not a nitrogen.";
        throw gromos::Exception(argv[0], msg.str());
      }
      // get the bonded neighbors
      Neighbours neighbours(sys, atoms.mol(i), atoms.atom(i));
      Neighbours::const_iterator it = neighbours.begin(), to = neighbours.end();
      for (; it != to; ++it) {
        // check whether neighbor atom is a nitrogen
        if (sys.mol(atoms.mol(i)).topology().atom(*it).mass() == hydrogen_mass)
          hydrogens[i].insert(*it);
      }
      // every selected nitrogen should have a neighbouring hydrogen (or up to 3)
      if (hydrogens[i].empty() || hydrogens[i].size() > 3) {
        ostringstream msg;
        msg << "Atom " << atoms.toString(i) << " has no or too many "
                "neighbouring hydrogens.";
        throw gromos::Exception(argv[0], msg.str());
      }
    }

    // get simulation time
    Time time(args);

    // get the averaging window
    unsigned int winframe = args.getValue<unsigned int>("winframe", false, 1);

    // loop over all trajectories
    int num_frames = 0;
    vector<gmath::Stat<double> > S2_win(atoms.size()); // holds the window averages S2 values
    vector<gmath::Matrix> S2_sum(atoms.size()); // the running sum of the tensor
    vector<gmath::Matrix> S2_sumwin(atoms.size()); // " window

    ofstream ts("OPts.out");
    ofstream tsw("OPwints.out");

    // setup the time series
    tsw.setf(ios::floatfield, ios::fixed);
    tsw.setf(ios::right, ios::adjustfield);
    tsw.precision(4);
    // print gromos numbers to title
    ts << "#" << setw(14) << "time";
    tsw << "#" << setw(14) << "time";
    for (unsigned int i = 0; i < atoms.size(); ++i) {
      ts << setw(10) << (atoms.gromosAtom(i) + 1);
      tsw << setw(10) << (atoms.gromosAtom(i) + 1);
    }
    ts << endl;
    tsw << endl;

    // define input coordinate
    for (Arguments::const_iterator iter = args.lower_bound("traj");
            iter != args.upper_bound("traj"); ++iter) {
      ic.open(iter->second);
      // loop over all frames
      while (!ic.eof()) {
        num_frames++;
        // load, gather and fit the system.
        ic >> sys >> time;
        (*pbc.*gathmethod)();
        if (rf != NULL)
          rf->fit(&sys);

        // write time to time series
        ts << time;

        // check whether we write out the window averaged values
        bool do_window = num_frames % winframe == 0;
        if (do_window)
          tsw << time;

        // calculate the z-vector between atom i-1 and i+1, normalize
        for (unsigned int i = 0; i < atoms.size(); i++) {
          // holds the N-H bond
          gmath::Vec nh; 
          set<int>::const_iterator it = hydrogens[i].begin();
          
          // in NMR it's not posisble to distinguish between the hydrogen bonds. 
          // so we calculate an average for secondary and tertiary NHx groups.
          switch (hydrogens[i].size()) {
            case 1:
            {
              // primary
              nh = (sys.mol(atoms.mol(i)).pos(*it) - atoms.pos(i)).normalize();
              break;
            }
            case 2:
            {
              // secondary
              int a = *it;
              int b = *(++it);
              nh = (((sys.mol(atoms.mol(i)).pos(a) - atoms.pos(i)).normalize() +
                      (sys.mol(atoms.mol(i)).pos(b) - atoms.pos(i)).normalize()) /
                      2.0).normalize();
              break;
            }
            case 3:
            {
              // tertiary
              int a = *it;
              int b = *(++it);
              int c = *(++it);
              nh = (((sys.mol(atoms.mol(i)).pos(a) - atoms.pos(i)).normalize() +
                      (sys.mol(atoms.mol(i)).pos(b) - atoms.pos(i)).normalize() +
                      (sys.mol(atoms.mol(i)).pos(c) - atoms.pos(i)).normalize()) /
                      3.0).normalize();
              break;
            }
            default:
              // just set N-H to zero. Could be conventient for prolines
              nh = gmath::Vec(0.0, 0.0, 0.0);
          }

          // calculate the sums for the averages
          gmath::Matrix dyadic_product(nh, nh);
          S2_sum[i] += dyadic_product;
          S2_sumwin[i] += dyadic_product;

          // write out the time series
          ts << setw(10) << S2_calc(S2_sum[i], num_frames);

          if (do_window) {
            const double s2 = S2_calc(S2_sumwin[i], winframe);
            tsw << setw(10) << s2;
            
            // save the value for statistics
            S2_win[i].addval(s2);

            // zero the matrix
            S2_sumwin[i] = gmath::Matrix();
          }
        } // for nitrogen atoms
        
        ts << endl;
        if (do_window) {
          tsw << endl;
        }
      } // while has frames
      ic.close();
    } // for trajectory
    
    // finish the time series
    ts.close();
    tsw.close();

    // print out results
    cout.setf(ios::right, ios::adjustfield);
    cout << "# MOL: molecule number" << endl
            << "# ATOM: atom number (in molecule)" << endl
            << "# ANAME: atom name" << endl
            << "# RES: residue number" << endl
            << "# RESNAME: residue name" << endl
            << "# S2: order parameter, averaged over whole trajectory" << endl
            << "# WINAV: order paramete, averaged over windows" << endl
            << "# WINRMSD: RMSD of WINAV" << endl
            << "# WINEE: error estimate of WINAV" << endl
            << "#" << endl
            << "# MOL" << setw(5) << "ATOM" << setw(8) << "ANAME" << setw(5) << "RES" 
            << setw(8) << "RESNAME" << setw(10) << "S2" << setw(10) << "WINAV"
            << setw(10) << "WINRMSD" << setw(10) << "WINEE" << endl;

    for (unsigned int i = 0; i < atoms.size(); ++i) {
      cout.setf(ios::floatfield, ios::fixed);
      cout.setf(ios::right, ios::adjustfield);
      cout.precision(4);
      cout << setw(5) << atoms.mol(i) + 1
              << setw(5) << atoms.atom(i) + 1
              << setw(8) << atoms.name(i)
              << setw(5) << atoms.resnum(i) + 1
              << setw(8) << atoms.resname(i)
              << setw(10) << S2_calc(S2_sum[i], num_frames)
              << setw(10) << S2_win[i].ave()
              << setw(10) << S2_win[i].rmsd()
              << setw(10) << S2_win[i].ee()
              << endl;
    }
    cout << endl;

    if (rf != NULL) {
      delete rf->getReference();
      delete rf;
    }

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

// calculate the double sum

double S2_calc(const gmath::Matrix & sum, int num_frames) {
  double double_sum = 0.0;
  for (unsigned int a = 0; a < 3; ++a) {
    for (unsigned int b = 0; b < 3; ++b) {
      const double ave = sum(a, b) / num_frames;
      double_sum += ave * ave;
    }
  }
  return 0.5 * (3.0 * double_sum - 1.0);
}
