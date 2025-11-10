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
 * @file make_sasa_top.cc
 * add sasa block to molecular topology file
 * required in order to use SASA/VOL implicit solvent model
 */

/**
 * @page programs Program Documentation
 *
 * @anchor make_sasa_top
 * @section make_sasa_top add sasa block to molecular topology file
 * @author @ref kb @ref ja
 * @date 23. 4. 2009
 *
 * Program make_sasa_top adds the atom-specific information required to use the
 * SASA/VOL implicit solvent model to the molecular topology file. It reads in
 * an existing molecular topology file created using @ref make_top, along with a
 * SASA/VOL specification library file, which contains the atom-specific SASA/VOL parameters.
 * The specification library file must be for the same force field as was used to create
 * the molecular topology file. The inclusion of hydrogen
 * atoms in the calculation of the sasa during the simulation can also be specified.
 * 
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@sasaspec</td><td>&lt;sasa specification library file&gt; </td></tr>
 * <tr><td> [\@noH</td><td>&lt;do not include hydrogen atoms (default: include)&gt;] </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
   make_sasa_top
     @topo       ex.top
     @sasaspec   sasaspec45b3.lib
   @endverbatim

 * <hr>
 */

#include <cassert>
#include <cstdlib>
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/LinearTopology.h"
#include "../src/gio/OutTopology.h"
#include "../src/gio/Ginstream.h"
#include "../src/gromos/Exception.h"

using namespace gcore;
using namespace gio;
using namespace args;

using namespace std;

struct sasa_parameter {
  double radius;
  double probability;
  double sigma;
};

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "sasaspec" << "noH";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo       <molecular topology file>\n";
  usage += "\t@sasaspec   <sasa specification file>\n";
  usage += "\t[@noH       <do not include hydrogen atoms (default: include)\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    gcore::System sys(it.system());

    // check whether we want hydrogens or not
    bool noH = false;
    if (args.count("noH") != -1)
      noH = true;

    // check if we have a sasaspec file
    if (args.count("sasaspec") != 1)
      throw gromos::Exception("sasaspec", "No sasa specification file");

    map<int, sasa_parameter> sasa_spec;
    {
      Ginstream spec_file(args["sasaspec"]);
      vector<string> buffer;
      spec_file.getblock(buffer);
      if (buffer[0] != "SASASPEC")
        throw gromos::Exception("sasaspec",
          "sasa specification file does not contain a SASASPEC block!");
      if (buffer[buffer.size() - 1].find("END") != 0)
        throw gromos::Exception("sasaspec", "sasa specification file " + spec_file.name() +
          " is corrupted. No END in SASASPEC"
          " block. Got\n"
          + buffer[buffer.size() - 1]);

      vector<string>::const_iterator it = buffer.begin() + 1, to = buffer.end() - 1;

      for (; it != to; ++it) {
        sasa_parameter param;
        istringstream line(*it);
        unsigned int num;
        line >> param.radius >> param.probability >> param.sigma >> num;
        if (line.fail())
          throw gromos::Exception("sasaspec",
            "bad line in SASASPEC block!");

        for (unsigned int i = 0; i < num; ++i) {
          int iac;
          line >> iac;
          if (line.fail()) {
            ostringstream msg;
            msg << "bad line in SASASPEC block: could not read " << num
                << " IACs from line.";
            throw gromos::Exception("sasaspec", msg.str());
          }
          sasa_spec[iac - 1] = param;
        } // for iacs
      } // SASASPEC block
    }

    // get total number of solute atoms to consider for SASA
    unsigned int numSASAatoms = 0;
    for (int m = 0; m < sys.numMolecules(); ++m) {
      if (!noH) {
        // if we're including H's, we can just add all atoms for this molecule
        numSASAatoms += sys.mol(m).numAtoms();
      } else {
        // check if each atom is an H or not
        for (int i = 0; i < sys.mol(m).numAtoms(); ++i) {
          if (!sys.mol(m).topology().atom(i).isH()) {
            numSASAatoms += 1;
          }
        }
      }
    }

    // write out original topology
    OutTopology ot(cout);
    ot.setTitle(it.title());
    ot.write(sys, it.forceField());

    // write out sasa block
    cout << "SASAPARAMETERS" << endl;
    cout << "#NRSASAA\n";
    cout << numSASAatoms << endl;
    cout << "#ISASA    RADI      PI     SIGMAI" << endl;

    for (int m = 0; m < sys.numMolecules(); ++m) {
      for (int i = 0; i < sys.mol(m).numAtoms(); ++i) {

        // if want hydrogens or this atom isn't a hydrogen
        if (!noH || !sys.mol(m).topology().atom(i).isH()) {

          // get iac
          int myiac = sys.mol(m).topology().atom(i).iac();
          map<int, sasa_parameter>::const_iterator result = sasa_spec.find(myiac);
          if (result == sasa_spec.end()) {
            ostringstream out;
            out << "No SASA parameters for atom type: " << myiac + 1;
            throw gromos::Exception("sasaspec", out.str());
          }

          const sasa_parameter & s = result->second;
          cout.precision(3);
          cout << setw(6) << i + 1;
          cout << setw(8) << s.radius << setw(8) << s.probability <<
          setw(11) << s.sigma << endl;
        }
      }
    }
    cout << "END" << endl;
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
