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

// gio_InG96.t.cc
#include "InG96.h"

#include <cassert>
#include <cstdlib>
#include <string>
#include <iostream>

#include "InTopology.h"
#include "OutG96S.h"
#include "../gcore/System.h"

using namespace std;
using namespace gcore;
using namespace gio;

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
  cout << "read topology and created system" << endl;

  InG96 ic;
  ic.open(file);
  ic.select("ALL");
  OutG96S oc;

  oc.open(cout);
  oc.select("ALL");
  oc.writeTitle(ic.title());

  //while(!ic.eof()){
  ic >> sys;
  cout << "read in something " << endl;

  oc << sys;
  //}
  return 0;
}


