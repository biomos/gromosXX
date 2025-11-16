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
#include "OutGromacs.h"

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <string>

#include "InTopology.h"
#include "InG96.h"
#include "../gromos/Exception.h"
#include "../gcore/GromosForceField.h"
#include "../gcore/System.h"

using namespace std;

using namespace gcore;
using namespace gio;

int main(int argc, char *argv[]) {
  if (argc != 2 && argc != 3) {
    cerr << "Usage: " + string(argv[0]) + " <Topology>   [<coordinate file>]\n";
    exit(1);
  }
  try {
    InTopology it(argv[1]);
    System sys(it.system());
    InG96 ic;
    if (argc == 3) {
      ic.open(argv[2]);
      ic.select("ALL");
      ic >> sys;
    }

    GromosForceField gff(it.forceField());
    OutGromacs ot(cout);
    ot.setTitle(it.title());
    ot.write(sys, gff);

    return 0;
  }  catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}

