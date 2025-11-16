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
 * @file grid_dens.cc
 * Calculate densities of atoms on a grid.
 */
/**
 * @page Contrib Program Documentation
 *
 * @anchor grid_dens
 * @section grid_dens Calculate densities of atoms on a grid.
 * @author @ref ns
 * @date 30.11.2010
 *
 * Program grid_dens is used to compute a density grid for the occurance of
 * atoms. The dimensions of the grid have to be specified along the three axes.
 * The specified atoms' positions are put in the box and are mapped on the grid
 * using the current box of the frame. So for NPT trajectories the densities are
 * to be interpreted with caution.
 * The grid is written to the standard output in Gaussian's Cube format which
 * can be read in VMD. Because Cube files have to contain atoms, a dummy atom
 * is placed at the origin.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> \@dim</td><td>&lt;dimensions of the grid&gt; </td></tr>
 * <tr><td> \@atoms</td><td>&lt;atoms to consider&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
 @topo  ex.top
 @pbc   r
 @dim   50 50 50
 @atoms 1:a
 @traj  ex.trc
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "../src/gmath/Vec.h"
#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InG96.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/debug.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Mesh.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;
using namespace gmath;

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "dim" << "atoms" << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc   <boundary type>\n";
  usage += "\t@dim   <dimensions of the grid, default 100 100 100>\n";
  usage += "\t@atoms <atoms for density>\n";
  usage += "\t@traj  <trajectory files>\n";

  try {
    Arguments args(argc, argv, knowns, usage);
    InTopology intopo(args["topo"]);
    System sys(intopo.system());

    // declare existence of reference system
    System refSys(intopo.system());
    // parse boundary conditions for refSys
    Boundary *pbc = BoundaryParser::boundary(refSys, args);

    // get the atoms. Check after frame was read.
    AtomSpecifier atoms(sys);
    for (Arguments::const_iterator it = args.lower_bound("atoms"),
            to = args.upper_bound("atoms"); it != to; ++it) {
      atoms.addSpecifier(it->second);
    }

    // get the mesh dimensions
    std::vector<unsigned int> dim = args.getValues<unsigned int>("dim", 3, false,
            Arguments::Default<unsigned int>() << 100 << 100 << 100);

    // initalize the mesh
    Mesh<double> grid(dim[0], dim[1], dim[2]);
    grid = 0.0;

    // use this to normalize the density
    double sum = 0.0;

    // loop over frames
    InG96 ic;
    for (Arguments::const_iterator
      iter = args.lower_bound("traj"),
            to = args.upper_bound("traj");
            iter != to; ++iter) {

      // open file
      ic.open((iter->second).c_str());
      ic.select("ALL");

      // loop over single trajectory
      while (!ic.eof()) {
        ic >> sys;
        if (!sys.hasBox || sys.box().ntb() == Box::truncoct) {
          throw gromos::Exception(argv[0], "Only works for rectangular or triclinic boxes!");
        }
        if (atoms.empty())
          throw gromos::Exception(argv[0], "No atoms given!");

        // use the current box for mapping. NPT?
        grid.setBox(sys.box());
        // get the centre of the box.
        Vec centre = (sys.box().K() + sys.box().L() + sys.box().M()) * 0.5;
        for(int i = 0; i < atoms.size(); ++i) {
          // put the particle on the positive box.
          Vec pos = pbc->nearestImage(centre, atoms.pos(i), sys.box());
          // increase occurance.
          grid(pos) += 1.0;
          sum += 1.0;
        }
      }
    }

    // normalize the grid, but avoid nans.
    if (sum != 0.0) {
      for (int i = 0; i < grid.numPoints(); ++i)
        grid(i) /= sum;
    }
    // write it out
    grid.write(cout);
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;


}
