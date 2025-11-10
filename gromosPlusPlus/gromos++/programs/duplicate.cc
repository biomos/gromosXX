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
 * @file duplicate.cc
 * Perform a search for duplicated molecules/atoms (and optionally
 * write out system without duplicated coordinates)
 */

/**
 * @page programs Program Documentation
 *
 * @anchor duplicate
 * @section duplicate Perform a search for duplicated molecules/atoms
 * @author @ref ns @ref ff
 * @date 23-3-07
 * 
 * Program duplicate searches for duplicated atoms, i.e. atoms having the
 * same coordinates as another atoms. If requested, a coordinate file with the
 * duplicated molecules removed is written out.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;input topology file&gt; </td></tr>
 * <tr><td> \@pos</td><td>&lt;input coordinate file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;periodic boundary conditions&gt; </td></tr>
 * <tr><td>[\@write</td><td>&lt;write out duplicate-filtered coordinates&gt;]</td></tr>
 * </table>
 *
 *
 * Example using a specification file:
 * @verbatim
  duplicate
    @topo    ref.topo
    @pos     exref.g96
    @pbc     r
    @write
 @endverbatim
 * <hr>
 */

#include <cassert>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <iostream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/bound/Boundary.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Vec.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace gmath;
using namespace bound;
using namespace args;

int main(int argc, char **argv) {
  Argument_List knowns;
  knowns << "topo" << "pos" << "pbc" << "write";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo       <input topology file>\n";
  usage += "\t@pos         <input coordinate file>\n";
  usage += "\t@pbc         <periodic boudary conditions>\n";
  usage += "\t[@write      <write out duplicate-filtered coordinate file, flag>]\n";

  // prepare cout for formatted output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(8);

  try {
    Arguments args(argc, argv, knowns, usage);
    // Original System
    // read topology
    args.check("topo", 1);
    InTopology it(args["topo"]);
    System sys1(it.system());
    // Red in file
    InG96 ic;
    ic.open(args["pos"]);
    ic.select("ALL");
    ic >> sys1;
    ic.close();

    // copy system (for double-loop)
    System sys2(sys1);

    // create empty system/outstream for output (optional)
    System outsys;
    outsys.box()=sys1.box();
    OutG96S oc;
    ostringstream title;
    // write titles
    // Parse boundary conditions
    Boundary & pbc = *BoundaryParser::boundary(sys1, args);
    bool write = args.count("write") != -1;
    vector<int> delmol;
    if (write) {
      title << "Duplicated atoms (mol1:atom1 <-> mol2:atom2)" << "\n";
      title << "-------------------------------------------" << "\n";
    } else {
      cout << "\nDuplicated atoms (mol1:atom1 <-> mol2:atom2)" << "\n";
      cout << "--------------------------------------------" << "\n";
    }
    // double loop (every atom with every atom) SOLUTE
    const Vec nullvec(0.0, 0.0, 0.0);

    for (int i = 0; i < sys1.numMolecules(); i++) {
      bool addmol = true;
      for (int j = 0; j < sys1.mol(i).numAtoms(); j++) {
        for (int k = i; k < sys2.numMolecules(); k++) {
          for (int l = j; l < sys1.mol(k).numAtoms(); l++) {
            // check if exactly the same coordinates (->duplicate)
            if (sys1.mol(i).pos(j) - pbc.nearestImage(sys1.mol(i).pos(j),
                    sys2.mol(k).pos(l), sys1.box()) == nullvec && (j != l || i != k)) {
              // write duplicates
              addmol = false;
              if (write) {
                title << "solute: " << i+1 << ":" << j+1 << "  <->  " << k+1 << ":" << l+1 << "\n";
              } else {
                cout << "solute: " << i+1 << ":" << j+1 << "  <->  " << k+1 << ":" << l+1 << "\n";
              }
            }
          }
        }
      }
      // optional molecule copy if no duplicated atom in there
      if (addmol) {
        if (write) {
          outsys.addMolecule(sys1.mol(i));
        }
      } else {
        delmol.push_back(i+1);
      }
    }
    // Run for Solvent
    int numdelsol = 0;
    for (int i = 0; i < sys1.numSolvents(); i++) {
      int solsize = sys1.sol(i).topology().numAtoms();
      // add empty solvent to outsys (fill later)
      const Solvent sol(sys1.sol(i).topology());
      outsys.addSolvent(sol);
      for (int j = 0; j < sys1.sol(i).numAtoms(); j += solsize) {
        bool addsol=true;
        for (int k = j; k < j + solsize; k++) {
          for (int l = i; l < sys2.numSolvents(); l++) {
            for (int m = j; m < sys2.sol(l).numAtoms(); m += solsize) {
              for (int n = m; n < m + solsize; n++) {
                if (sys1.sol(i).pos(k) - pbc.nearestImage(sys1.sol(i).pos(k),
                        sys2.sol(l).pos(n), sys1.box()) == nullvec && (k != n)) {
                  // write duplicates
                  addsol = false;
                  if (write) {
                    title << "solvent: " << j / solsize << ":" << k << "  <->  " << m / solsize << ":" << n << "\n";
                  } else {
                    cout << "solvent: " << j / solsize << ":" << k << "  <->  " << m / solsize << ":" << n << "\n";
                  }
                }
              }
            }
          }
        }
        // optional solvent-molecule copy if no duplicated atom in there
        if (addsol) {
          if (write) {
            for (int o = j; o < j + solsize; o++) {
              outsys.sol(i).addPos(sys1.sol(i).pos(o));
            }
          }
        } else {
          numdelsol++;
        }
      }
    }

    // output of filtered system (optional)
    if (write) {
      title << "\nDeleted Solute-Molecules:  " << setw(6) << delmol.size() << "\n";
      for(vector<int>::const_iterator it = delmol.begin(), to = delmol.end(); it != to; ++it)
        title << "  - Molecule " << *it << "\n";
      title << "Deleted Solvent-Molecules: " << setw(6) << numdelsol;
      oc.open(cout);
      oc.select("ALL");
      oc.writeTitle(title.str());
      oc << outsys;
    } else {
      cout << "\nDuplicated Solute-Molecules:  " << setw(6) << delmol.size() << "\n";
      for(vector<int>::const_iterator it = delmol.begin(), to = delmol.end(); it != to; ++it)
        cout << "  - Molecule " << *it << "\n";
      cout << "Duplicated Solvent-Molecules: " << setw(6) << numdelsol << "\n";
    }
    return 0;

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
}
