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
 * @file nscatt.cc
 */

/**
 * @page programs Program Documentation
 *
 * @anchor nscatt
 * @section nscatt short description
 * @author A. Eichenberger
 * @date March 2011
 *
 * Description...
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * ...
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
 nscatt
    @topo   ex.top
    @traj   ex.tr
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <string>

#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/System.h"
#include "../src/utils/NeutronScattering.h"
#include "../src/gromos/Exception.h"


using namespace std;
using namespace args;
using namespace gcore;
using namespace gio;
using namespace utils;

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "atoms" << "cut" << "grid" << "scattlen" << "sigma" << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t@centre   <AtomSpecifier: atoms to be considered as centre atoms>\n";
  usage += "\t@with     <AtomSpecifier: atoms to be considered as with atoms>\n";
  usage += "\t[@grid    <number of data points>]\n";


  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology
    args.check("topo", 1);
    InTopology it(args["topo"]);
    System sys(it.system());

    args.check("scattlen", 1);

    if(args.count("traj") < 1) {
      throw gromos::Exception(argv[0], "no trajectory file(s) specified");
    }
    NS ns(&sys, &args);

    // set the centre and with atoms
    if(args.count("atoms") < 1) {
      throw gromos::Exception(argv[0], "no atoms defined (@atoms)");
    }
    {
      Arguments::const_iterator iter = args.lower_bound("atoms");
      Arguments::const_iterator to = args.upper_bound("atoms");
      for (; iter != to; iter++) {
        ns.addAtoms(iter->second.c_str());
      }
    }

    // this is the sequence it HAS to be done:
    //   1) get the number of combinations (also resetting all the vector lengths
    //   2) set the system (for all subvectors too
    //   3) set the atoms to the AtomSpecifiers (for all subvectors)
    ns.getCombinations();
    ns.setSystem(&sys);
    ns.setRDFatoms();
    ns.getWeights();

    // set the grid number to the specified integer number, if there is an @ grid flag
    if(args.count("grid") > 0) {
      stringstream ss;
      ss << args["grid"];
      int grid;
      ss >> grid;
      if(ss.fail() || ss.bad()) {
        stringstream msg;
        msg << "could not convert " << args["grid"] << " into an integer number";
        throw gromos::Exception(argv[0], msg.str());
      }
      ns.setGrid(grid);
    }

    if(args.count("cut") > 0) {
      stringstream ss;
      ss << args["cut"];
      double cut;
      ss >> cut;
      if(ss.fail() || ss.bad()) {
        stringstream msg;
        msg << "could not convert " << args["cut"] << " into an double";
        throw gromos::Exception(argv[0], msg.str());
      }
      ns.setCut(cut);
    }

    ns.readScattlen(args["scattlen"]);
    ns.readSigma(args["sigma"]);
    ns.check();

    ns.print(cerr);
    ns.calcRDFsInterAll();
    ns.calcSintra();
    ns.calcSinter();
    ns.printSintra(cout);
    ns.printSinter(cout);
    ns.calcIntensity();
    ns.printRDFs(cout);
    ofstream fout("intensity.dat");
    ns.printIntensity(fout);
    ns.printS(cerr);
    fout.close();

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

