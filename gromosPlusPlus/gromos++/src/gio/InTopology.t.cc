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
#include "InTopology.h"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>

#include "../gromos/Exception.h"
#include "OutTopology.h"
#include "../gcore/GromosForceField.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"

using namespace gcore;
using namespace gio;

using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " + std::string(argv[0]) + " <Topology>\n";
    exit(1);
  }
  try {
    cout << "create an it" << endl;

    InTopology it(argv[1]);

    cout << "done with it" << endl;
    cout << "create a system" << endl;

    System sys(it.system());
    cout << "done with sys" << endl;
    cout << sys.mol(0).topology().numRes() << endl;

    GromosForceField gff(it.forceField());
    std::cout << sys.numMolecules() << endl;
    OutTopology ot(std::cout);
    ot.setTitle(it.title());
    ot.write(sys, gff);

    return 0;
  }  catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}
