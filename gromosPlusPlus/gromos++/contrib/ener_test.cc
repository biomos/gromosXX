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
 * @file ener_test.cc
 * This is a program to test the Energy.cc class which will act as a template
 * class for a force class (still to be implemented).
 * 
 * @page contrib Contrib program documentation
 *
 * @anchor ener_test
 * @section calculates the energy between two groups of
 *          atoms specified by atom specifiers
 * @author @ref ae
 * @date November 12, 2012
 *
 * Program ener_test ist a temporary program only. It is used to test the
 * Energy.cc class by calculating the interaction energy of two groups of atoms
 * specified by atom specifiers. The result will be compared to the energy
 * trajectory of a short simulations to make sure the two energies as calculated
 * from MD++ and GROMOS++ are the same.
 * 
 * If the energy class seems to be fine, it will be used as a template for a
 * force class to be implemented in GROMOS++.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;periodic boundary conditions [and gather option]&gt; </td></tr>
 * <tr><td> \@atomsA</td><td>&lt;@ref AtomSpecifier atoms of group A (atom specifier)&gt; </td></tr>
 * <tr><td> \@atomsB</td><td>&lt;@ref AtomSpecifier atoms of group B (atom specifier)&gt; </td></tr>
 * <tr><td> \@trc</td><td>&lt;positional simulation trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  rdf
    @topo    ex.top
    @pbc     r cog
    @atomsA  1:1-10
    @atomsB  1:11-12
    @trc     simulation.trc.gz
 @endverbatim
 *
 * <hr>
 */

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <iomanip>
#include <string>

#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/System.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/Energy.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/bound/Boundary.h"
#include "../src/args/GatherParser.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/utils/groTime.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace args;
using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "atoms" << "trc";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t@pbc      <periodic boundary conditions [and gather option]>\n";
  usage += "\t@atoms    <atoms of group (atom specifier)>\n";
  //usage += "\t@atomsA   <atoms of group A (atom specifier)>\n";
  //usage += "\t@atomsB   <atoms of group B (atom specifier)>\n";
  usage += "\t@trc      <positional simulation trajectory files>\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology
    args.check("topo", 1);
    InTopology it(args["topo"]);
    System sys(it.system());
    System refSys(it.system()); // a reference topology, in our case the same as the actual topology

    // The GROMOS force field as read from the topology
    GromosForceField gff(it.forceField());

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    //parse gather method
    //Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, refSys, args);

    // parameters of the non-bonded interactions
    double cutOff = 1.4;
    double eps = 61;
    double kappa = 0;

    // create Time object
    Time time(args);

    // add the atoms of groups A and B to the atom specifiers
    AtomSpecifier atoms(sys);
    //AtomSpecifier atomsA(sys);
    //AtomSpecifier atomsB(sys);
    {
      Arguments::const_iterator start = args.lower_bound("atoms");
      Arguments::const_iterator end = args.upper_bound("atoms");
      for (; start != end; start++) {
        atoms.addSpecifier(start->second);
      }
      /*
      Arguments::const_iterator iter = args.lower_bound("atomsA");
      Arguments::const_iterator to = args.upper_bound("atomsA");
      for (; iter != to; iter++) {
        atomsA.addSpecifier(iter->second.c_str());
      }
      iter = args.lower_bound("atomsB");
      to = args.upper_bound("atomsB");
      for (; iter != to; iter++) {
        atomsB.addSpecifier(iter->second);
      }
       */
    }


    // the energy class to be used
    Energy en(sys, gff, *pbc);

    // the standard parameters to be overwritten if specified in the input file

    en.setAtoms(atoms); // select all the (solute) atoms of the system
    en.setCutOff(cutOff);
    en.setPairlistType("CHARGEGROUP");
    en.setRF(eps, kappa);
    en.setRFexclusions(true); // as we decided to always do (remember the 
    // GROMOSCOMPAT block issue)

    // just some output to check the program
    cerr << "# number of atoms selected: " << atoms.size() << endl;
    //cerr << "# number of atoms in atom group A: " << atomsA.size() << endl;
    //cerr << "# number of atoms in atom group B: " << atomsB.size() << endl;
    
    cout << setw(15) << "time" << setw(15) << "total"
            << setw(15) << "non-bonded"
            << setw(15) << "Coul./RF"
            << setw(15) << "vdW\n";
    
    // loop over the trajectory files
    for (Arguments::const_iterator iter = args.lower_bound("trc");
            iter != args.upper_bound("trc"); ++iter) {

      // define input coordinates
      InG96 ic;
      ic.open(iter->second);
      ic.select("ALL");

      // loop over all frames
      while (!ic.eof()) {
        ic >> sys >> time;
        // calculate the energies
        en.calc();
        cout.precision(9);
        cout << setw(15) << time.time() << setw(15) << en.tot()
                << setw(15) << en.nb()
                << setw(15) << en.el()
                << setw(15) << en.vdw()
                << endl;
      }
    }

    /*

    // read in cut-off distance
    Rdf.setCut(args.getValue<double>("cut", true));

    // read in grid number
    Rdf.setGrid(args.getValue<int>("grid", true));

    // Check if intramolecular rdf is included
    bool nointra = false;
    if (args.count("nointra") >= 0) nointra = true;

    // calculate the rdf
    if(nointra) {
      Rdf.calculateInter();
    } else {
      Rdf.calculateAll();
    }

    // print the result
    Rdf.print(cout);
     */
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

