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
#include "CheckTopo.h"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>

#include "../gio/InTopology.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gio/InG96.h"
#include "../gromos/Exception.h"

using namespace gcore;
using namespace gio;
using namespace utils;

using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "Usage: " + string(argv[0]) + " <Topology>\n";
    exit(1);
  }
  try {
    InTopology it(argv[1]);
    System sys(it.system());

    CheckTopo ct(sys.mol(0).topology());
    ct.checkBonds();
    ct.checkAngles();
    ct.checkImpropers();
    ct.checkChargeGroups();
    ct.checkExclusions();

    cout << ct.numErrors() << endl;
    for (int i = 0; i < ct.numErrors(); i++) {
      cout << ct.error(i) << endl;
    }
    return 0;
  }  catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}
