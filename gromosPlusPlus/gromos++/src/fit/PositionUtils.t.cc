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

// fit_PositionUtils.t.cc

#include "PositionUtils.h"

#include <cassert>
#include <cstdlib>
#include <string>
#include <iostream>

#include "../gio/InG96.h"
#include "../gcore/System.h"
#include "../gio/InTopology.h"
#include "../gio/OutG96.h"
#include "../gmath/Vec.h"

using namespace gcore;
using namespace gio;
using namespace fit;
using namespace gmath;
using namespace std;

int debug_level = 0;

int main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << "Usage: " + string(argv[0]) + " <Topology> <G96inates>\n";
    exit(1);
  }
  string top = argv[1];
  string file = argv[2];
  // read Simulation data

  InTopology it(top);
  System sys(it.system());

  InG96 ic;
  ic.open(file);


  OutG96 oc;
  oc.open(cout);
  oc.writeTitle(ic.title());

  while (!ic.eof()) {
    ic >> sys;
    Vec com = PositionUtils::com(sys);
    cout << "Center of Mass: " << com[0] << ' ' << com[1]
            << ' ' << com[2] << endl;
    Vec cog = PositionUtils::cog(sys);
    cout << "Center of Geometry: " << cog[0] << ' ' << cog[1]
            << ' ' << cog[2] << endl;
    PositionUtils::translate(&sys, -com);
    cout << "Mol(0) fitted to center of mass:\n";
    oc << sys;

  }
  return 0;
}
