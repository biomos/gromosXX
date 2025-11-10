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

// utils_Neighbours.t.cc
#include "Neighbours.h"

#include <cassert>
#include <cstdlib>
#include <string>
#include <iostream>

#include "../gcore/System.h"
#include "../gio/InTopology.h"

using namespace gcore;
using namespace gio;
using namespace utils;

using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 4) {
    cerr << "Usage: " << argv[0] << " <topology> <mol nr.> <atom nr.>" << endl;
    exit(1);
  }
  string top = argv[1];
  int mol = atoi(argv[2]) - 1;
  int nr = atoi(argv[3]) - 1;

  InTopology it(top);
  System sys = it.system();

  Neighbours neigh(sys, mol, nr);
  cout << "Neighbours of atom " << nr + 1 << " of molecule " << mol + 1 << ": ";
  for (Neighbours::iterator iter = neigh.begin(); iter != neigh.end(); ++iter)
    cout << *iter + 1 << ' ';
  cout << endl;
  return 0;
}
