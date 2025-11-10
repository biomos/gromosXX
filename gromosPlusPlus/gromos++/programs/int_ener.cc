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
 * @file int_ener.cc
 * Recalculates interaction energies between two specified sets of atoms
 */

/**
 * @page programs Program Documentation
 *
 * @anchor int_ener
 * @section int_ener Recalculate interaction energies between two specified sets of atoms
 * @author @ref ja
 * @date 09-07-2009
 *
 * Program int_ener recalculates the nonbonded interaction energy between
 * two non-overlapping sets of solute atoms (A and B) using the interaction parameters
 * specified in the molecular topology file. If there is overlap between A and B,
 * the overlapping atoms will automatically be excluded. It can also compute the interaction
 * energy between a specified group of solute atoms (A) and the solvent. If
 * a time-series is requested, the total nonbonded interaction is printed at each
 * time point, along with the van der Waals and electrostatic contributions.
 *
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> \@atomsA</td><td>&lt;first group of @ref AtomSpecifier "atoms"&gt; </td></tr>
 * <tr><td> [\@atomsB</td><td>&lt;second group of @ref AtomSpecifier "atoms"&gt;] </td></tr>
 * <tr><td> [\@solvent</td><td>&lt;compute energy between atomsA and solvent&gt;] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> [\@timeseries</td><td>&lt;print time-series&gt;] </td></tr>
 * <tr><td> [\@timespec</td><td>&lt;timepoints at which to compute the energy: ALL (default), EVERY or SPEC (if time-series)&gt;] </td></tr>
 * <tr><td> [\@timepts</td><td>&lt;timepoints at which to compute the energy (if time-series and timespec EVERY or SPEC)&gt;] </td></tr>
 * <tr><td> [\@cut</td><td>&lt;cut-off distance (default: 1.4)&gt;] </td></tr>
 * <tr><td> [\@eps</td><td>&lt;epsilon for reaction field contribution (default: 1.0)&gt;] </td></tr>
 * <tr><td> [\@kap</td><td>&lt;kappa for reaction field contribution (default: 0.0)&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;position trajectory file(s)&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  int_ener
    @topo       ex.top
    @pbc        r
    @atomsA     1:3-13
    @atomsB     1:20-30
    @time       0 0.2
    @timeseries
    @cut        1.4
    @eps        61.0
    @kap        0.0
    @traj       ex.trj
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/Energy.h"
#include "../src/utils/groTime.h"
#include "../src/gromos/Exception.h"


using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

// function for skipping frames
bool computeEne(int i, std::string const & timespec, vector<int> const & timepts,
        unsigned int & timesWritten, bool & done);

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "atomsA" << "atomsB" << "solvent" << "time" << "timeseries" <<
          "timespec" << "timepts" << "cut" << "eps" << "kap" << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo         <molecular topology file>\n";
  usage += "\t@pbc          <boundary type> [<gather method>]\n";
  usage += "\t@atomsA       <first group of atoms>\n";
  usage += "\t[@atomsB       <second group of atoms>]\n";
  usage += "\t[@solvent      <compute energy between atomsA and solvent>]\n";
  usage += "\t[@time        <time and dt>]\n";
  usage += "\t[@timeseries  <print time-series>]\n";
  usage += "\t[@timespec    <timepoints at which to compute the energy: ALL (default), EVERY or SPEC (if time-series)>]\n";
  usage += "\t[@timepts     <timepoints at which to compute the energy (if time-series and timespec EVERY or SPEC)>]\n";
  usage += "\t[@cut         <cut-off distance (default: 1.4)>]\n";
  usage += "\t[@eps         <epsilon for reaction field correction (default: 1.0)>]\n";
  usage += "\t[@kap         <kappa for reaction field correction (default: 0.0)>]\n";
  usage += "\t@traj         <position trajectory file(s)>\n";


  try {
    Arguments args(argc, argv, knowns, usage);

    //  read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    GromosForceField gff(it.forceField());

    System refSys(it.system());

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args); 

    // declare the energy class
    Energy en(sys, gff, *pbc);

    //  set first group of atoms and start total group
    AtomSpecifier atomsA(sys);
    AtomSpecifier atomsAB(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("atomsA");
      Arguments::const_iterator to = args.upper_bound("atomsA");
      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        atomsA.addSpecifier(spec);
        atomsAB.addSpecifier(spec);
      }
    }

    // check whether the interaction is with solvent or with atomsB
    bool solv_ene = false;
    AtomSpecifier atomsB(sys);
    if (args.count("solvent") >= 0) {
      solv_ene = true;
      // just set up energy for atomsA
      en.setAtoms(atomsA);
      if (args.count("atomsB") > 0) {
        throw gromos::Exception("int_ener",
                "Cannot use @solvent and @atomsB simultaneously.\n");
      }
    } else if (args.count("atomsB") > 0) {

      // set second group of atoms
      {
        Arguments::const_iterator iter = args.lower_bound("atomsB");
        Arguments::const_iterator to = args.upper_bound("atomsB");
        for (; iter != to; iter++) {
          string spec = iter->second.c_str();
          atomsB.addSpecifier(spec);
          atomsAB.addSpecifier(spec);
        }
      }
      // set up all atoms for energy calculation (energy can only handle one AtomSpecifier)
      en.setAtoms(atomsAB);
    } else {
      throw gromos::Exception("int_ener",
              "Must specify either @solvent or @atomsB.\n");
    }

    // set non-bonded parameters
    //   get cut-off distance
    {
      double cutoff = args.getValue<double>("cut", false, 1.4);
      en.setCutOff(cutoff);
    }
    
    //  get epsilon and kappa
    {
      double eps = args.getValue<double>("eps", false, 1.0);
      double kap = args.getValue<double>("kap", false, 0.0);
      en.setRF(eps, kap);
    }

    // get simulation time if given
    Time time(args);

    // parse timespec
    string timespec = "ALL";
    vector<int> timepts;
    if (args.count("timespec") > 0) {
      timespec = args["timespec"];
      if (timespec != "ALL" && timespec != "EVERY" && timespec != "SPEC")
        throw gromos::Exception("int_ener",
              "timespec format " + timespec + " unknown.\n");
      if (timespec == "EVERY" || timespec == "SPEC") {
        for (Arguments::const_iterator it = args.lower_bound("timepts");
                it != args.upper_bound("timepts"); ++it) {
          int bla = atoi(it->second.c_str());
          timepts.push_back(bla);
        }
        if (timepts.size() == 0) {
          throw gromos::Exception("int_ener",
                  "if you give EVERY or SPEC you have to use "
                  "@timepts as well");
        }
        if (timepts.size() != 1 && timespec == "EVERY") {
          throw gromos::Exception("int_ener",
                  "if you give EVERY you have to give exactly"
                  " one number with @timepts");
        }
      }
    }

    // check if we want a time-series of the energies
    bool ene_ts = false;
    if (args.count("timeseries") >= 0) {
      ene_ts = true;

      // print title of time-series
      cout << "# Time"
              << "                   vdw"
              << "                    elec"
              << "                 Total"
              << endl;
    }

    // declare some variables for averaging
    double ave_vdw = 0.0;
    double ave_elec = 0.0;
    double ave_tot = 0.0;

    // set counters for timespecs
    int numTimepoints = 0;
    // number of frames that have been written
    unsigned int timesWritten = 0;
    // for SPEC: so that we stop trying when all requested timepoints are written
    bool done = false;

    // define input coordinate
    InG96 ic;

    // loop over all trajectories
    for (Arguments::const_iterator
      iter = args.lower_bound("traj"),
            to = args.upper_bound("traj");
            iter != to; ++iter) {

      // open file
      ic.open((iter->second).c_str());
      ic.select("ALL");

      // loop over single trajectory
      while (!ic.eof()) {

        ic >> sys >> time;
        (*pbc.*gathmethod)();

        // update number for timespec
        numTimepoints++;

        // check whether to use this block or move on
        if (computeEne(numTimepoints, timespec, timepts, timesWritten, done)) {

          // declare variables for this frame
          double vdw = 0.0;
          double elec = 0.0;


          // if we want energy with solvent:
          if (solv_ene) {
            // calculate the energies (for atomsA)
            en.calc();
            // loop over selected atoms
            for (unsigned int i = 0; i < atomsA.size(); ++i) {
              vdw += en.vdw_s(i);
              elec += en.el_s(i);
            }
          }
          // otherwise, do pairwise energy for atoms A and B
          else {

            // loop over first set of selected atoms (A)
            for (unsigned int i = 0; i < atomsA.size(); ++i) {

	      // loop over the combined set of atoms (AB)
              for (unsigned int j = atomsA.size(); j < atomsAB.size(); ++j) {

		// set totals to zero
                double vdw_ij = 0.0;
                double elec_ij = 0.0;

                // check for overlapping selections (by molecule and atom numbers)
                if (!(atomsA.mol(i) == atomsAB.mol(j) &&
                        atomsA.atom(i) == atomsAB.atom(j))) {

                  // compute pair energy for atoms i and jj
                  en.calcPair(i, j, vdw_ij, elec_ij);

		// overlapping atom selections: not an error, simply that they are excluded
                //} else {
		//	// why wasn't this printed?
		//	cout << "Overlapping atoms " << atomsA.atom(i) << " " << atomsAB.atom(j) << endl;
                }

                // add atom pair energies to sums for this frame
                vdw += vdw_ij;
                elec += elec_ij;

              } // end loop j atomsAB
            } // end loop i atomsA
          } // end if pairwise A:B

          // total interaction of selected atoms with solvent
          double tot = vdw + elec;

          // print ouput if time-series
          if (ene_ts) {
            cout.precision(10);
            cout.setf(ios::right, ios::adjustfield);
            cout << setw(6) << time
                    << setw(22) << vdw
                    << setw(22) << elec
                    << setw(22) << tot
                    << endl;
          }

          //store some averages
          ave_vdw += vdw;
          ave_elec += elec;
          ave_tot += tot;

        } // end if computeJval
        if (done)
          break;
      } // ic loop
      ic.close();
    } // end loop over trajectories

    // print out averages
    if (timesWritten > 1) {
      cout.precision(10);
      cout.setf(ios::right, ios::adjustfield);
      cout << endl << "# Averages               vdw"
              << "                    elec                 Total" << endl << "#"
              << setw(38) << ave_vdw / timesWritten
              << setw(22) << ave_elec / timesWritten
              << setw(22) << ave_tot / timesWritten
              << endl;
    }

    // case where there is only one structure
    else if (timesWritten == 1) {
      cout.precision(10);
      cout.setf(ios::right, ios::adjustfield);
      cout << endl << "# Energies               vdw"
              << "                    elec                 Total" << endl << "#"
              << setw(38) << ave_vdw / timesWritten
              << setw(22) << ave_elec / timesWritten
              << setw(22) << ave_tot / timesWritten
              << endl;
    }


  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

bool computeEne(int i, std::string const & timespec, vector<int> const & timepts,
        unsigned int & timesWritten, bool & done) {
  if (timespec == "ALL") {
    ++timesWritten;
    return true;
  } else if (timespec == "EVERY" && i % timepts[0] == 0) {
    ++timesWritten;
    return true;
  } else if (timespec == "SPEC") {
    for (unsigned int j = 0; j < timepts.size(); ++j) {
      if (timepts[j] == i) {
        ++timesWritten;
        if (timesWritten == timepts.size())
          done = true;
        return true;
      } // compute ene?
    } // times
  }
  return false;
}
