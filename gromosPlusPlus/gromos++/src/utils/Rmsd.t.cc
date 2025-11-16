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

// utils_Rmsd.t.cc
#include "Rmsd.h"

#include <cassert>
#include <cstdlib>
#include <string>
#include <iostream>

#include "../fit/RotationalFit.h"
#include "../fit/Reference.h"
#include "../gio/InG96.h"
#include "../gcore/System.h"
#include "../gio/InTopology.h"

using namespace gcore;
using namespace gio;
using namespace utils;
using namespace fit;

using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 4) {
    cerr << "Usage: " + string(argv[0]) + " <Topology> <RefG96inates> <G96inates>\n";
    exit(1);
  }
  string top = argv[1];
  string refFile = argv[2];
  string file = argv[3];
  // read Simulation data

  InTopology it(top);

  System refsys(it.system());

  InG96 ic(refFile);
  ic >> refsys;
  ic.close();

  System sys(refsys);

  Reference ref(&refsys);

  ref.addClass(0, "CA");
  ref.addClass(0, "N");
  ref.addClass(0, "C");

  RotationalFit rf(&ref);
  Rmsd rmsd(&ref);

  ic.open(file);

  while (!ic.eof()) {
    ic >> sys;
    rf.fit(&sys);
    double d = rmsd.rmsd(sys);
    cout << d << endl;
  }
  ic.close();
  return 0;
}
