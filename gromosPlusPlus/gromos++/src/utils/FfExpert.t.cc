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
#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>

#include "FfExpert.h"
#include "../gio/InBuildingBlock.h"
#include "../gromos/Exception.h"

using namespace gcore;
using namespace gio;
using namespace utils;

using namespace std;

int debug_level = 0;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "Usage: " + string(argv[0]) + " <BuildingBlock>\n";
    exit(1);
  }
  try {
    InBuildingBlock ibb(argv[1]);

    FfExpert exp(ibb.building());

    cout << "succes" << endl;
    vector<FfExpert::counter> v;
    exp.name2iac("C", v);

    cout << "we have " << v.size() << " atoms with C" << endl;
    for (unsigned int i = 0; i < v.size(); i++) {
      cout << v[i].type << " " << v[i].occurence << endl;
    }
    exp.iac2mass(12, v);

    cout << "we have " << v.size() << " masses for iac=12" << endl;
    for (unsigned int i = 0; i < v.size(); i++) {
      cout << v[i].type << " " << v[i].occurence << endl;
    }
    int b = 2;

    exp.iac2charge(b, v);
    cout << "we have " << v.size() << " charges for iac=" << b << endl;
    for (unsigned int i = 0; i < v.size(); i++) {
      cout << exp.charge(v[i].type) << " " << v[i].occurence << endl;
    }


  }  catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}
