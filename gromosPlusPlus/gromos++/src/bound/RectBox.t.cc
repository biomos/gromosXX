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

// bound_RectBox.t.cc

#include "RectBox.h"

#include <cassert>
#include <cstdlib>
#include <string>
#include <iostream>

#include "../gio/InG96.h"
#include "../gcore/System.h"
#include "../gio/InTopology.h"
#include "../gio/OutG96.h"
#include "../gcore/System.h"
#include "../gcore/Solvent.h"

using namespace gcore;
using namespace gio;
using bound::RectBox;

using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << "Usage: " + string(argv[0]) + " <Topology> <Filename>\n";
    exit(1);
  }
  string top = argv[1];
  string file = argv[2];

  // read Simulation data

  InTopology it(top);
  System sys(it.system());
  Solvent sol = sys.sol(0);

  InG96 ic;


  ic.open(file);
  ic.select("ALL");
  cout << "sys.numSolvents: " << sys.numSolvents() << endl;
  cout << "sol.numCoords: " << sol.numPos() << endl;
  OutG96 oc;

  ic >> sys;
  cout << "sys.numSolvents: " << sys.numSolvents() << endl;
  RectBox pbc(&sys);
  cout << "sys.numSolvents after RectBox pbc(&sys): " << sys.numSolvents() << endl;

  oc.open(cout);
  oc.select("ALL");
  oc.writeTitle(ic.title());

  while (!ic.eof()) {
    pbc.coggather();
    oc << sys;
    ic >> sys;
  }
  pbc.coggather();
  oc << sys;
  return 0;
}
