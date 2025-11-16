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
#include "InParameter.h"

#include <cassert>
#include <iostream>
#include <cstdlib>
#include <string>

#include "../gromos/Exception.h"
#include "OutTopology.h"
#include "../gcore/GromosForceField.h"
#include "../gcore/System.h"

using namespace std;
using namespace gcore;
using namespace gio;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "Usage: " + string(argv[0]) + " <Parameterfile>\n";
    exit(1);
  }
  try {
    InParameter ip(argv[1]);
    GromosForceField gff(ip.forceField());
    System sys;

    OutTopology ot(cout);
    ot.setTitle(ip.title());
    ot.write(sys, gff);

    return 0;
  }  catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}
