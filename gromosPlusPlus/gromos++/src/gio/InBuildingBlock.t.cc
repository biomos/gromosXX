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
#include "InBuildingBlock.h"

#include <cassert>
#include <iostream>
#include <cstdlib>
#include <string>

#include "../gcore/BuildingBlock.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/BbSolute.h"
#include "../gcore/SolventTopology.h"
#include "../gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "Usage: " + string(argv[0]) + " <Building Block file>\n";
    exit(1);
  }
  try {
    InBuildingBlock ibb(argv[1]);
    BuildingBlock bld(ibb.building());

    cout << "Number of BbSolutes: " << bld.numBbSolutes() << endl;
    if (bld.numBbSolutes() > 23) {

      cout << " BbSolute number 24: " << bld.bb(23).resName() << endl;
      cout << "  has " << bld.bb(23).numAtoms() << " atoms and "
              << bld.bb(23).numPexcl() << " preceding exclusions" << endl;
    }
    cout << "Number of BbSolvents: " << bld.numBbSolvents() << endl;
    if (bld.numBbSolvents() > 3) {

      cout << "  BbSolvent number 3: " << bld.bs(2).solvName() << endl;
      cout << "  has " << bld.bs(2).numAtoms() << " atoms" << endl;
    }

    cout << "Number of BbEnds: " << bld.numBbEnds() << endl;
    if (bld.numBbEnds() > 2) {

      cout << "  BbEnd number 2: " << bld.be(2).resName() << endl;
      cout << "  will replace " << bld.be(2).rep() << " atoms" << endl;
    }

    int index = bld.findBb("DADE");
    DihedralIterator di(bld.bb(index - 1));
    int count = 0;
    for (; di; ++di) count++;
    cout << "DADE has " << count << " dihedrals" << endl;
    return 0;
  }  catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}





