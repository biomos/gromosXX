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
 * @file rdf.cc
 * calculates a radial distribution function
 */

/**
 * @page programs Program Documentation
 *
 * @anchor rdf
 * @section rdf calculates a radial distribution function
 * @author @ref co
 * @date 28.7.2006
 *
 * Program rdf calculates radial distribution functions over structure files or
 * trajectories. The radial distribution function, g(r), is defined here as the
 * probability of finding a particle of type J at distance r from a central
 * particle I relative to the same probability for a homogeneous distribution
 * of particles J around I. Program rdf calculates g(r) for a number of
 * discreet distances r(k), separated by distance dr as
 *
 * @f[ g(r) = \frac{N_J(k)}{4\pi r^2 dr \rho_J} @f]
 *
 * where @f$N_J(k)@f$ is the number of particles of type J found at a distance 
 * between r(k) - 1/2 dr and r(k) + 1/2 dr and @f$\rho_J@f$ is the number 
 * density of particles J. If particles I and J are of the same type, 
 * @f$\rho_J@f$ is corrected for that. At long distances, g(r) will generally 
 * tend to 1. 
 *
 * Both atoms of type I and J can be solute atoms, solvent atoms as well as 
 * @ref VirtualAtom "virtual atoms". If more than one particle of type I is
 * specified, rdf calculates the average radial distribution function for all
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
 * <tr><td> [\@nointra</td><td>&lt;skip all intramolecular contributions&gt;] </td></tr>
 * <tr><td> [\@doDCF</td><td>&lt;calculate the dipole-dipole correlations&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  rdf
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
#include "../src/gromos/Exception.h"
#include "../src/utils/RDF.h"


using namespace std;
using namespace args;
using namespace gcore;
using namespace gio;
using namespace utils;


int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "centre" << "with"
          << "cut" << "grid" << "nointra" << "traj" << "pbc" << "doDCF";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t@pbc        <boundary type> [<gather method>]\n";
  usage += "\t@centre <atoms to take as centre>\n";
  usage += "\t@with   <atoms to calculate distances for>\n";
  usage += "\t@cut    <maximum distance>\n";
  usage += "\t@grid   <number of points>\n";
  usage += "\t[@nointra   <skip intramolecular atoms>]\n";
  usage += "\t[@doDCF     <compute dipole-dipole correlation> <norm>]\n";
  usage += "\t@traj   <trajectory files>";


  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology
    args.check("topo", 1);
    InTopology it(args["topo"]);
    System sys(it.system());

    // the class for the calculation of the rdf
    RDF Rdf(&sys, &args);

    // set the centre and with atoms
    {
      Arguments::const_iterator iter = args.lower_bound("centre");
      Arguments::const_iterator to = args.upper_bound("centre");
      for (; iter != to; iter++) {
        Rdf.addCenters(iter->second.c_str());
      }
      iter = args.lower_bound("with");
      to = args.upper_bound("with");
      for (; iter != to; iter++) {
        Rdf.addWiths(iter->second.c_str());
      }
    }

    // read in cut-off distance
    Rdf.setCut(args.getValue<double>("cut", true));

    // read in grid number
    Rdf.setGrid(args.getValue<int>("grid", true));

    // check if we want to compute the dipole-dipole correlation function
    if(args.count("doDCF") >=0 ) {
      Rdf.setDCF(true);
      Rdf.DCFnorm(false);
      if(args.count("doDCF") >= 1 && args["doDCF"] == "norm"){
        // cerr << "Normalizing the DCF" << endl;
        Rdf.DCFnorm(true);
      }
    }

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

    if(args.count("doDCF") >= 0)
	Rdf.print_DCF(cout);

  }  catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

