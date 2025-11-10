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
 * @file rmsd_cmp.cc
 * calculates atom-positional root-mean-square deviations
 */
/**
 * @page programs Program Documentation
 *
 * @anchor rmsd_cmp
 * @section rmsd_cmp calculates atom-positional root-mean-square deviations between two trajectories
 * @author @ref nb
 * @date 26.7.06
 *
 * The structural deformation of a molecule with respect to a reference
 * structure can be expressed in terms of a root-mean-square deviation (rmsd)
 * of the position of selected atoms. Program rmsd_cmp calculates the rmsd over two
 * molecular trajectory. If requested a least-square rotational fit is performed before
 * to the rmsd calculation. The fit
 * can be performed using a different set of atoms than the calculation of the 
 * rmsd. If no fit is required "no" should be given.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gathermethod&gt;] </td></tr>
 * <tr><td> \@time</td><td>&lt;@ref utils::Time "time and dt"&gt; </td></tr>
 * <tr><td> \@atomsrmsd</td><td>&lt;@ref AtomSpecifier "atoms" to consider for rmsd&gt; </td></tr>
 * <tr><td> [\@atomsfit</td><td>&lt;@ref AtomSpecifier "atoms" to consider for fit&gt;] </td></tr>
 * <tr><td> \@reftraj</td><td>&lt;reference trajectory &gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  ditrans
    @topo       ex.top
    @pbc        r
    @time       0 0.1
    @atomsrmsd  1:CA
    @atomsfit   1:CA,C,N
    @reftraj    exref.tr
    @traj       ex.tr

    @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/args/ReferenceParser.h"
#include "../src/utils/Rmsd.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/gromos/Exception.h"


using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;
using namespace std;
using namespace gmath;

int main(int argc, char **argv){
  Argument_List knowns; 
  knowns << "topo" << "traj" << "atomsfit" << "atomsrmsd" << "prop" << "pbc" << "reftraj"
         << "time" << "debug" << "fit";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo       <molecular topology file>\n";
  usage += "\t@pbc        <boundary type> [<gathermethod>]\n";
  usage += "\t@time       <time and dt>\n";
  usage += "\t@atomsrmsd  <atoms to consider for rmsd>\n";
  usage += "\t[@atomsfit  <atoms to consider for fit>]\n";
  usage += "\t[@prop      <properties>\n";
  usage += "\t@reftraj    <reference coordinates (if absent, the first frame of @traj is reference)>\n";
  usage += "\t@traj       <trajectory files>\n";


  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);

  try {
    Arguments args(argc, argv, knowns, usage);

    // get simulation time
    Time time(args);
    // read topology
    InTopology it(args["topo"]);
    System refSys(it.system());
    //System sys(refSys);
    System sys(it.system());

    // Parse boundary conditions
    Boundary *refpbc = BoundaryParser::boundary(refSys, args);
    // GatherParser
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, refSys, args);

    if (args.count("reftraj") <= 0)
      throw gromos::Exception("rmsd_cmp", "No '@reftraj' argument given!");

    //cout << " ref mol " << refSys.numMolecules() << " solv " << refSys.sol(0).numPos() << endl;

    int debug = 0;
    if (args.count("debug") > 0)
      debug = 1;
    int fit = 0;
    if (args.count("fit") > 0)
      fit = atoi(args.lower_bound("fit")->second.c_str());

    // System for calculation
    Reference refrmsd(&refSys);
    AtomSpecifier fitatoms(refSys);
    AtomSpecifier rmsdatoms(sys);

    //get rmsd atoms
    if (args.count("atomsrmsd") > 0) {
      {
        Arguments::const_iterator iter = args.lower_bound("atomsrmsd");
        Arguments::const_iterator to = args.upper_bound("atomsrmsd");

        for (; iter != to; iter++) {
          string spec = iter->second.c_str();
          rmsdatoms.addSpecifier(spec);
        }
      }
      if (rmsdatoms.size() == 0)
        throw gromos::Exception("rmsd", "No rmsd-atoms specified!");

      refrmsd.addAtomSpecifier(rmsdatoms);
    }

    //try for fit atoms
    if (args.count("atomsfit") > 0) {
      Arguments::const_iterator iter = args.lower_bound("atomsfit");
      Arguments::const_iterator to = args.upper_bound("atomsfit");

      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        fitatoms.addSpecifier(spec);
      }
    } else {
      cout << "# @atomsrmsd atoms are taken for fit." << endl;
      const vector<string> & spec = rmsdatoms.toString();

      for (vector<string>::const_iterator it = spec.begin(), to = spec.end();
              it != to; ++it) {
        fitatoms.addSpecifier(*it);
      }
    }

    // Parse boundary conditions for sys
    Boundary *pbc = BoundaryParser::boundary(sys, args);

    // Property container also for reference system
    PropertyContainer prop_ref(refSys, pbc);
    PropertyContainer prop_sys(sys, pbc);
    {
      std::string prop;
      Arguments::const_iterator iter = args.lower_bound("prop");
      Arguments::const_iterator to = args.upper_bound("prop");
      // we read in all properties specified by the user

      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        prop += " " + spec;
      }

      prop_ref.addSpecifier(prop);
      prop_sys.addSpecifier(prop);
    }

    if (args.count("atomsrmsd") < 0 && args.count("prop") < 0) {
      throw gromos::Exception("rmsd", "No rmsd atoms or property specified!");
    }

    //Vec cog=PositionUtils::cog(refSys, reffit);

    RotationalFit * rf = NULL;
    if (fitatoms.size()) {
      Reference * ref = new Reference(&refSys);
      ref->addAtomSpecifier(fitatoms);
      rf = new RotationalFit(ref);
    }

    Rmsd rmsd(&refrmsd);
    rmsd.addproperty(&prop_ref, &prop_sys);

    int numFrames = 0;

    //cout << rmsd.rmsd(refSys) << endl;

    InG96 ic;
    InG96 icref;

    // loop over all trajectories
    for (Arguments::const_iterator iter = args.lower_bound("traj"), 
            iter_ref = args.lower_bound("reftraj");
            iter != args.upper_bound("traj"); ++iter, ++iter_ref) {
      ic.open(iter->second);
      icref.open(iter_ref->second);

      // loop over all frames
      while (!ic.eof()) {

        numFrames++;
        
        // Normal Trajectory
        ic.select("ALL");
        ic >> sys >> time;
        if (!sys.hasPos)
          throw gromos::Exception("rmsd",
                "Unable to read POSITION(RED) block from "
                "trajectory file.");
        (*pbc.*gathmethod)();

        // Reference Trajectory
        icref.select("ALL");
        icref >> refSys;
        if (!refSys.hasPos)
          throw gromos::Exception("rmsd",
                "Unable to read POSITION(RED) block from "
                "reference positions file.");
        (*refpbc.*gathmethod)();
    
        if (fitatoms.size())
          rf->fit(&sys);

        double r = rmsd.rmsd(sys);
        double rprop = rmsd.rmsdproperty(sys);

        //	cout.precision(2);
        //	cout << time;
        //	cout.precision(5);
        if (args.count("atomsrmsd") > 0 && args.count("prop") > 0) {
          cout.precision(2);
          cout << time;
          cout.precision(9);
          cout << setw(15) << r << setw(15) << rprop << endl;
        } else if (args.count("atomsrmsd") > 0) {
          cout.precision(2);
          cout << time;
          cout.precision(5);
          cout << setw(10) << r << endl;
        } else if (args.count("prop") > 0) {
          cout.precision(2);
          cout << time;
          cout.precision(9);
          cout << setw(15) << rprop << endl;
        }
      } // for every frame
      ic.close();
      icref.close();
    } // for every trajectory

    if (rf != NULL) {
      delete rf->getReference();
      delete rf;
    }

  }  catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


