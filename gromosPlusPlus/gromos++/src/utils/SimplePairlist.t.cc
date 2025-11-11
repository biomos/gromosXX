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

#include <cassert>
#include <set>
#include <iostream>
#include <cstdlib>
#include "../gio/InTopology.h"
#include "AtomSpecifier.h"
#include "SimplePairlist.h"
#include "../gcore/System.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Molecule.h"
#include "../gcore/Solvent.h"
#include "../bound/Boundary.h"
#include "../bound/Vacuum.h"
#include "../gio/InG96.h"

using namespace gcore;
using namespace gio;
using namespace utils;

using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 4) {
    cerr << "Usage: " + string(argv[0]) + " <Topology> <atomspecifier> <coordinates>\n";
    exit(1);
  }
  try {
    InTopology it(argv[1]);
    System sys(it.system());
    string s = argv[2];

    AtomSpecifier as(sys, s);
    bound::Boundary *pbc = new bound::Vacuum(&sys);
    InG96 ic(argv[3]);
    ic.select("ALL");

    ic >> sys;
    ic.close();

    SimplePairlist pl(sys, *pbc, 0.8);
    pl.setAtom(as.mol(0), as.atom(0));
    pl.setType("ATOMIC");

    pl.calc();
    pl.removeExclusions();
    pl.remove14Exclusions();

    for (int i = 0; i < pl.size(); i++)
      cout << pl.mol(i) + 1 << ":" << pl.atom(i) + 1 << endl;



    return 0;
  } catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}
