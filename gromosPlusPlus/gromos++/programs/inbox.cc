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
 * @file inbox.cc
 * Put atoms into box according to pbc.
 */

/**
 * @page programs Program Documentation
 *
 * @anchor inbox
 * @section inbox put atoms into the box
 * @author @ref ns
 * @date 11-03-2009
 *
 * Program inbox puts the atoms into the positive computational box
 * according to the periodic boundary conditions. It can be used to visualize
 * the computational box in a crystal simulation. The connectivity and gathering
 * of chargegroups is ignored, thus the chargegroups (and solvent molecules)
 * won't be gathered after application of this program.
 *
 * Using the \@atoms argument one can specify the atoms which are put into the
 * box. All other atoms are not affected by the program. By default all atoms
 * are put into the box.
 *
 * Using the \@shift argument one can shift the computational box (with respect 
 * to the default positive box). The shifts in x, y and z direction are by 
 * default in nm but can also be chosen in terms of boxlengths, if keyword 
 * boxlength is given as last parameter.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;a trajectory&gt; </td></tr>
 * <tr><td>[\@atoms</td><td>&lt;@ref AtomSpecifier "atoms"&gt;]</td></tr>
 * <tr><td>[\@shift</td><td>&lt;shifts in x y z [boxlength]&gt;]</td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
 inbox
    @topo             ex.top
    @pbc              r
    @atoms            a:a
    @shift            -0.5 -0.5 -0.5 boxlength
    @traj             ex.tr
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <map>
#include <string>
#include <iostream>
#include <sstream>

#include "../src/gcore/Box.h"
#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/args/OutformatParser.h"
#include "../src/gromos/Exception.h"

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace gmath;
using namespace utils;
using namespace std;

int main(int argc, char **argv) {
  Argument_List knowns;
  knowns << "topo" << "pbc" << "atoms" << "traj" << "shift" << "outformat";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo           <molecular topology file>\n";
  usage += "\t@pbc            <boundary type> \n";
  usage += "\t[@shift         <x y z <boxlength> shift in nm (default) or boxlengths>\n";
  usage += "\t[@atoms         <atoms>]\n";
  usage += "\t[@outformat <output coordinates format>]\n";
  usage += "\t@traj           <trajectory files>\n";


  try {
    Arguments args(argc, argv, knowns, usage);

    InTopology it(args["topo"]);
    System sys(it.system());
    
    string ext;
    OutCoordinates *oc = OutformatParser::parse(args, ext);

    Boundary *pbc = BoundaryParser::boundary(sys, args);

    AtomSpecifier atoms(sys);
    {
      Arguments::const_iterator to = args.upper_bound("atoms"),
              iter = args.lower_bound("atoms");
      if (iter == to) {
        // add all atoms
        atoms.addSpecifier("a:a");
        atoms.addSpecifier("s:a");
      } else {
        for (; iter != to; iter++)
          atoms.addSpecifier(iter->second);
      }
    }
    
    bool shift=false;
    bool boxlength=false;
    gmath::Vec shifts(0.,0.,0.);
    if (args.count("shift")>=0) {
      std::string argstrings="";
      std::string bl;
      shift=true;
      Arguments::const_iterator iter=args.lower_bound("shift");
      Arguments::const_iterator to=args.upper_bound("shift");
      unsigned int n=0;
      while (iter != to) {
        argstrings+=" ";
        argstrings+=iter->second;
        n+=1;
        iter++;
      }
      std::istringstream is(argstrings);
      if (n>=3)
        is >> shifts[0] >> shifts[1] >> shifts[2];
      if (n>=4) {
        is >> bl;
        if (bl == "boxlength") boxlength = true;
        else cerr << "\nWARNING: optional fourth parameter for @shift can\n"
                  << "         only be 'boxlength', otherwise defaulting to nm!!";
      }  
      if (n>4 || n<3)
        throw gromos::Exception("inbox",
				"invalid number of arguments for @shift!");
    }

    // loop over all trajectories
    InG96 ic;
    Arguments::const_iterator iter = args.lower_bound("traj"),
            to = args.upper_bound("traj");

    if (iter == to)
      throw gromos::Exception(argv[0], "no trajetctory given.");

    // prepare the output
    oc->open(cout);
    oc->writeTitle("Trajectory put in the box");

    for (;iter != to; ++iter) {
      // open file
      ic.open((iter->second).c_str());

      // loop over single trajectory
      while (!ic.eof()) {
        ic.select("ALL");
        ic >> sys;

        if (!sys.hasBox)
          throw gromos::Exception(argv[0], "no box in frame.");

        // get the centre of the box
        const Vec centre = (sys.box().K() * 0.5) + (sys.box().L() * 0.5) +
                (sys.box().M() * 0.5);

        for(unsigned int i = 0; i < atoms.size(); ++i) {
          if (shift) {
            if (boxlength) {
              gmath::Vec boxshift=(sys.box().K() * shifts[0]) + (sys.box().L() * shifts[1]) +
                (sys.box().M() * shifts[2]);
              atoms.pos(i) = atoms.pos(i)  + boxshift;
            }
            else {
              atoms.pos(i) = atoms.pos(i)  + shifts;
            }
          }
            atoms.pos(i) = atoms.pos(i) - pbc->nearestImage(atoms.pos(i), centre, sys.box()) + centre;
        }

        // write out the frame
        oc->select("ALL");
        *oc << sys;
      } // while frames
      ic.close();
    } // for traj

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

