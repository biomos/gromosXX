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
 * @file local_comp.cc
 * calculates a radial distribution function
 */

/**
 * @page programs Program Documentation
 *
 * @anchor local_comp
 * @section local_comp calculates a radial distribution function
 * @author @ref co
 * @date 28.7.2006
 *
 * Program local_comp calculates local composition of a binary mixture over structure files or
 * trajectories. The local fractions, @f[x_{ii}(r)@f], is defined as 
 * 
 * @f[ x_{ii}(r) = \frac{n_{ii}(r)}{n_{ji}(r) + n_{ii}(r)} @f].
 * 
 * The number of atoms of type i or j around atoms of type i is obtained from
 * 
 * @f[ n_{ij}(r) = \int_0^r \frac{N_i}{V} g_{ij}(r') 4\pi r^2 dr' @f]
 * 
 * Program local_comp calculates @f[x_{ii}(r)@f] for a number of
 * discreet distances r(k), separated by distance dr. 
 *
 * Both atoms of type i and j can be solute atoms, solvent atoms as well as 
 * @ref VirtualAtom "virtual atoms". If more than one particle of type i is
 * specified, local_comp calculates the average radial distribution function for all
 * specified atoms.
 *
 * The boundary type of the read systems is read from the GENBOX block of the
 * trajectory files.
 *
 * This program is parallelised.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@centre</td><td>&lt;@ref AtomSpecifier "atoms" to take as centre&gt; </td></tr>
 * <tr><td> \@with</td><td>&lt;@ref AtomSpecifier "atoms" to calculate distances for&gt; </td></tr>
 * <tr><td> \@cut</td><td>&lt;maximum distance&gt; </td></tr>
 * <tr><td> \@grid</td><td>&lt;number of points&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  local_comp
    @topo   ex.top
    @pbc    r
    @centre 1:45
    @with   s:OW
    @cut    3.0
    @grid   100
    @traj   ex.tr
 @endverbatim
 *
 * <hr>
 */

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <string>

#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/System.h"
#include "../src/utils/RDF.h"
#include "../src/gromos/Exception.h"


using namespace std;
using namespace args;
using namespace gcore;
using namespace gio;
using namespace utils;


int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "centre" << "with"
          << "cut" << "grid" << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t@centre <atoms to take as centre>\n";
  usage += "\t@with   <atoms to calculate distances for>\n";
  usage += "\t@cut    <maximum distance>\n";
  usage += "\t@grid   <number of points>\n";
  usage += "\t@traj   <trajectory files>\n";


  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology
    args.check("topo", 1);
    InTopology it(args["topo"]);
    System sys(it.system());

    // first and last trajectory
    Arguments::const_iterator firsttrj = args.lower_bound("traj");
    Arguments::const_iterator lasttrj = args.upper_bound("traj");

    // the class for the calculation of the rdf
    RDF rdf(&sys, &args);

    // set the center and with atoms
    {
      Arguments::const_iterator iter = args.lower_bound("centre");
      Arguments::const_iterator to = args.upper_bound("centre");
      for (; iter != to; iter++) {
        rdf.addCenters(iter->second.c_str());
      }
      iter = args.lower_bound("with");
      to = args.upper_bound("with");
      for (; iter != to; iter++) {
        rdf.addWiths(iter->second.c_str());
      }
    }

    // read in cut-off distance
    rdf.setCut(args.getValue<double>("cut", true));

    // read in grid number
    rdf.setGrid(args.getValue<int>("grid", true));

    rdf.calculateLocal();

    // print the result
    rdf.printLocal(cout);

  }  catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

