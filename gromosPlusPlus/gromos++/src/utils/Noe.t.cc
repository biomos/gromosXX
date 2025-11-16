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

// utils_Noe.t.cc
#include "Noe.h"

#include <cassert>
#include <cstdlib>
#include <string>
#include <iostream>

#include "../gcore/System.h"
#include "../gio/InTopology.h"
#include "../gio/InG96.h"
#include "../gromos/Exception.h"

using namespace utils;
using namespace gcore;
using namespace gio;

using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 4) {
    cerr << "Usage: " << argv[0] << " <topology> <coordinates> <Noe line>" << endl;
    exit(1);
  }
  string top = argv[1];
  string coord = argv[2];

  string line = argv[3];

  try {

    InTopology it(top);
    System sys = it.system();

    InG96 ic(coord);
    ic >> sys;


    Noe noe(sys, line, 0.1, 0.15);

    cout << "The distance restraints are \n";
    for (int i = 0; i < noe.numDistances(); ++i)
      cout << noe.distRes(i) << endl;

    cout << "\nThe distances calculated are:\n";
    for (int i = 0; i < noe.numDistances(); ++i)
      cout << noe.distance(i) << ' ';

    cout << "\n\nThe following reference distances were found. \nUncorrected: ";
    for (int i = 0; i < noe.numReferences(); ++i)
      cout << noe.reference(i) << ' ';
    cout << "\nCorrected: ";
    for (int i = 0; i < noe.numReferences(); ++i)
      cout << noe.correctedReference(i) << ' ';
    cout << endl;


  } catch (gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }

  return 0;
}
