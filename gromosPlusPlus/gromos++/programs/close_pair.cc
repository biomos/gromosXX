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
 * @file close_pair.cc
 * find atom pairs of two molecules that are close to each other
 */

/**
 * @page programs Program Documentation
 *
 * @anchor close_pair
 * @section close_pair
 * @author @ref dw
 * @date 01. 07. 2010
 *
 * Program close_pair finds the closest atom pairs of two molecules.
 * The searching for molecule i loops over all (or specified) solute atoms of the other molecules,
 * and if the distance between the mol(i).atom(j) and mol(k).atom(l) is shorter
 * than a value (dc default 0.3 nm),
 * the searching for the atom pairs between molecules i - k is stopped. Periodicity is
 * also considered if the pbc is defined.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@groupA</td><td>&lt;@ref AtomSpecifier "atoms" of molecules to be gathered&gt; </td></tr>
 * <tr><td> \@groupB</td><td>&lt;@ref AtomSpecifier "atoms" of molecules to be as reference&gt; </td></tr>
 * <tr><td> \[\@dist</td><td>&lt;lower limit of distance. default 0.3 nm]&gt; </td></tr>
 * <tr><td> \[\@time</td><td>&lt;t0 and dt]&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
 close_pair
    @topo       ex.top
    @pbc        r
    @groupA   a:C
    @groupB   a:a
    [@dist     0.3  ]
    [@time     0  1 ]
    @traj       ex.tr
@endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <string>
#include <fstream>
#include <vector>
#include <iomanip>
#include <iostream>

#include <time.h>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Matrix.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/groTime.h"
#include "../src/gcore/Box.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gmath;
using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "traj" << "groupA" << "groupB" << "pbc" << "time" << "dist" << "debug";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo      <molecular topology file>\n";
  usage += "\t@pbc         <boundary type>\n";
  usage += "\t@groupA      <molecule group to be gathered>\n";
  usage += "\t@groupB      <molecule group to be reference for gathering of groupA>]\n";
  usage += "\t[@dist       <lower limit of distance. default 0.3 nm>]\n";
  usage += "\t[@time       <t0 and dt>]\n";
  usage += "\t@traj        <trajectory files>\n";


  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);


  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    InG96 ic;

    //get time
    utils::Time time(args);

    double dist = args.getValue<double>("dist", false, 0.3);
    // we will calc only the square of distance
    double dc = dist * dist;

    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);

    //get the atoms
    AtomSpecifier groupA(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("groupA");
      Arguments::const_iterator to = args.upper_bound("groupA");

      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        groupA.addSpecifier(spec);
      }
    }

    AtomSpecifier groupB(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("groupB");
      Arguments::const_iterator to = args.upper_bound("groupB");

      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        groupB.addSpecifier(spec);
      }
    }
    int debug = 0;
    if (args.count("debug") > 0) {
      debug = atoi(args["debug"].c_str());
    }
    if (debug == 1) {
      for (unsigned int i = 0; i < groupA.size(); ++i) {
        cout << "# groupA mol " << groupA.mol(i) << " size " << groupA.size() << endl;
      }
      for (unsigned int j = 0; j < groupB.size(); ++j) {
        cout << "# groupB mol " << groupB.mol(j) << " size " << groupB.size() << endl;
      }
    }

    if (groupA.size() == 0)
      throw gromos::Exception("groupA",
            "No atoms specified for groupA");

    groupA.sort();
    groupB.sort();

    double boxsquare = sys.box().X();
    boxsquare *= boxsquare / 4.0;

    cout << setw(10) << "#    Frame" << setw(8) << "molA"
            << setw(8) << "resA" << setw(6) << "atomA"
            << setw(8) << "molB" << setw(8) << "resB" << setw(6) << "atomB"
            << setw(10) << "atomA" << setw(8) << "atomB"
            << setw(12) << "dist [nm]" << endl;

    // Parse boundary conditions for sys
    pbc = BoundaryParser::boundary(sys, args);

    int numFrames = 0;

    //vector<Vec> apos;
    //Vec apos(0.,0.,0.);
    vector<double> apos2;
    vector<int> aresid, aatmid;
    vector<string> aatom;
    vector<int> bmolid, bresid, batmid;
    vector<string> batom;

    //dimension= sys.numMolecules() * sys.numMolecules();
    apos2.resize(sys.numMolecules(), 0.0); // shortest dist
    aresid.resize(sys.numMolecules(), 0); // molA:ResID
    aatmid.resize(sys.numMolecules(), 0); // molA:AtomID
    aatom.resize(sys.numMolecules(), ""); // molA:AtomName

    bmolid.resize(sys.numMolecules(), 0);
    bresid.resize(sys.numMolecules(), 0);
    batmid.resize(sys.numMolecules(), 0);
    batom.resize(sys.numMolecules(), "");

    // vector to save temporary values
    vector< vector<double> > buff_apos2(sys.numMolecules(), apos2);
    vector< vector<int> > buff_aresid(sys.numMolecules(), aresid), buff_aatmid(sys.numMolecules(), aatmid);
    vector< vector<string> > buff_aatom(sys.numMolecules(), aatom);
    vector< vector<int> > buff_bmolid(sys.numMolecules(), bmolid), buff_bresid(sys.numMolecules(), bresid), buff_batmid(sys.numMolecules(), batmid);
    vector< vector<string> > buff_batom(sys.numMolecules(), batom);

    double temp; // save the dist^2 and then compare with apos2, save the shorter one

    double clo;

    //loop over trajectory
    for (Arguments::const_iterator iter = args.lower_bound("traj");
            iter != args.upper_bound("traj"); ++iter) {
      ic.open(iter->second);
      // loop over all frames
      while (!ic.eof()) {
        // read frame
        ic.select("ALL");
        //ic >> sys >> time;
        ic >> sys;

        if (debug) {
          clo = clock();
        }

        int pass = 1000;

        for (unsigned int i = 0; i < groupA.size(); ++i) {
          // skip hydrogen atoms
          int m = groupA.mol(i);

          if (debug >= 3)
            cout << "# search for atom  " << i << "   " << m << ":" << groupA.resnum(i) << " " << groupA.name(i) << endl;

          if (groupA.resname(i) == "")
            throw gromos::Exception("groupA",
                  "atom name not specified properly.");
          if (m != pass && m > 0)
            if (groupA.mass(i) != 1.008) {
              for (unsigned int j = 0; j < groupB.size(); ++j) {
                if (debug >= 3) {
                  cout << "# groupB size " << groupB.size() << endl;
                  cout << "# molA " << groupA.mol(i)
                          << " atom " << i
                          << " name " << groupA.name(i)
                          << " mass " << groupA.mass(i)
                          << "\t molB " << groupB.mol(j)
                          << " atom " << j
                          << " name " << groupB.name(j)
                          << " mass " << groupB.mass(j) << endl;
                  cout << "# groupB.mol(j) " << groupB.mol(j) << " groupA.mol(i) " << groupA.mol(i) << endl;
                }
                if (groupB.name(j) == "")
                  throw gromos::Exception("groupB",
                        "atom name not specified properly.");
                if (groupB.mol(j) < groupA.mol(i)) {
                  if (groupB.mass(j) != 1.008) {
                    Vec ref = pbc->nearestImage(groupA.pos(i), groupB.pos(j), sys.box());
                    temp = (groupA.pos(i) - ref).abs2();
                    if (debug)
                      cout << "# done for nearest image " << j << endl;
                    if (debug) {
                      double distance = (groupA.pos(i) - groupB.pos(j)).abs2();
                      cout << m << " : " << i << "\t" << groupB.mol(j) << " : " << j
                              << "\t" << temp << "\t" << apos2[m] << "\t" << distance << endl;
                    }

                    if (debug >= 3)
                      cout << "# done for pair: " << m << ":" << i << " " << groupA.name(i) << "\t"
                      << groupB.mol(j) << ":" << j << " " << groupB.name(j) << "\t" << temp << endl;

                    if (buff_apos2[m][groupB.mol(j)] == 0 && temp != 0) {
                      buff_apos2[m][groupB.mol(j)] = temp;
                      buff_aresid[m][groupB.mol(j)] = groupA.resnum(i);
                      buff_aatmid[m][groupB.mol(j)] = groupA.atom(i);
                      buff_aatom[m][groupB.mol(j)] = groupA.name(i);

                      buff_bmolid[m][groupB.mol(j)] = groupB.mol(j);
                      buff_bresid[m][groupB.mol(j)] = groupB.resnum(j);
                      buff_batmid[m][groupB.mol(j)] = groupB.atom(j);
                      buff_batom[m][groupB.mol(j)] = groupB.name(j);
                      if (debug)
                        cout << "buff term m " << m << " molB " << groupB.mol(j) << " initialized " << endl;
                    } else {
                      if (buff_apos2[m][groupB.mol(j)] > temp) {
                        buff_apos2[m][groupB.mol(j)] = temp;
                        buff_aresid[m][groupB.mol(j)] = groupA.resnum(i);
                        buff_aatmid[m][groupB.mol(j)] = groupA.atom(i);
                        buff_aatom[m][groupB.mol(j)] = groupA.name(i);

                        buff_bmolid[m][groupB.mol(j)] = groupB.mol(j);
                        buff_bresid[m][groupB.mol(j)] = groupB.resnum(j);
                        buff_batmid[m][groupB.mol(j)] = groupB.atom(j);
                        buff_batom[m][groupB.mol(j)] = groupB.name(j);
                      }
                    }
                    if (apos2[m] == 0) {
                      apos2[m] = temp;
                      aresid[m] = groupA.resnum(i);
                      aatmid[m] = groupA.atom(i);
                      aatom[m] = groupA.name(i);

                      bmolid[m] = groupB.mol(j);
                      bresid[m] = groupB.resnum(j);
                      batmid[m] = groupB.atom(j);
                      batom[m] = groupB.name(j);

                      if (debug >= 3)
                        cout << " apos2 initialized : " << m << " : " << i << "\t" << groupB.mol(j) << " : " << j << "\t" << apos2[m] << endl;
                    } else {
                      if (apos2[m] > temp) {
                        if (debug >= 3)
                          cout << " redefine apos2 : " << m << " : " << i << "\t" << groupB.mol(j) << " : " << j << "\t" << temp << endl;
                        apos2[m] = temp;
                        aresid[m] = groupA.resnum(i);
                        aatmid[m] = groupA.atom(i);
                        aatom[m] = groupA.name(i);

                        bmolid[m] = groupB.mol(j);
                        bresid[m] = groupB.resnum(j);
                        batmid[m] = groupB.atom(j);
                        batom[m] = groupB.name(j);

                        if (apos2[m] < dc) {
                          pass = m;
                          break;
                        }
                      }
                    }
                  } else {
                    if (debug >= 3)
                      cout << "# this is an H atom, thus skipped. " << endl;
                  }
              }
              }
            }
        }
        for (unsigned int i = 0; i < buff_apos2.size(); ++i) {
          for (unsigned int j = 0; j < buff_apos2[i].size(); ++j)
            if (buff_apos2[i][j] != 0)
              cout << "# " << setw(8) << numFrames
                    << setw(5) << i + 1
                    << ":res(" << setw(5) << buff_aresid[i][j] + 1
                    << ":" << setw(5) << buff_aatom[i][j]
                    << ") " << setw(5) << buff_bmolid[i][j] + 1
                    << ":res(" << setw(5) << buff_bresid[i][j] + 1
                    << ":" << setw(5) << buff_batom[i][j]
                    << ") # " << setw(6) << buff_aatmid[i][j] + 1
                    << setw(8) << buff_batmid[i][j] + 1
                    << setw(12) << sqrt(buff_apos2[i][j]) << endl;
        }

        cout << endl << "# Summary" << endl << "# =======" << endl
                << setw(10) << "#    Frame" << setw(8) << "molA"
                << setw(8) << "resA" << setw(6) << "atomA"
                << setw(8) << "molB" << setw(8) << "resB" << setw(6) << "atomB"
                << setw(10) << "atomA" << setw(8) << "atomB"
                << setw(12) << "dist [nm]" << endl;
        for (unsigned int i = 0; i < apos2.size(); ++i) {
          if (apos2[i] != 0.0)
            cout << "# " << setw(8) << numFrames
                  << setw(5) << i + 1
                  << ":res(" << setw(5) << aresid[i] + 1
                  << ":" << setw(5) << aatom[i]
                  << ") " << setw(5) << bmolid[i] + 1
                  << ":res(" << setw(5) << bresid[i] + 1
                  << ":" << setw(5) << batom[i]
                  << ") # " << setw(6) << aatmid[i] + 1
                  << setw(8) << batmid[i] + 1
                  << setw(12) << sqrt(apos2[i]) << endl;
        }

        // generate a file containing list readable by program atominfo for lazy users
        ofstream atomlist("atominfo.atomspec");
        for (unsigned int i = 0; i < apos2.size(); ++i) {
          if (apos2[i] != 0.0)
            atomlist << setw(5) << i + 1
                  << ":res(" << aresid[i] + 1
                  << ":" << aatom[i]
                  << ") " << setw(5) << bmolid[i] + 1
                  << ":res(" << bresid[i] + 1
                  << ":" << batom[i]
                  << ") " << endl;
        }
        atomlist.close();


        numFrames++;
        if (debug) {
          cout << endl << "# time used (s) : " << (clock() - clo) / (CLOCKS_PER_SEC) << endl;
        }
      }
    }
    ic.close();
  }//end loop over trajectory
  catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

