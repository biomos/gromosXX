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

// gcore_Exclusion.t.cc

#include "Exclusion.h"

#include <iostream>

using gcore::Exclusion;
using namespace std;

static void output(const Exclusion &e) {
  cout << "size: " << e.size() << endl
          << "Exclusions: ";

  for (int i = 0; i < e.size(); ++i)
    cout << e.atom(i) << ' ';
  cout << endl;
}

int main() {
  Exclusion e;
  e.insert(3);
  e.insert(2);
  e.insert(5);
  e.insert(4);
  e.insert(4);

  output(e);

  Exclusion g = e;
  output(g);

  Exclusion f;
  f = g;
  output(f);
  f.erase(3);
  output(f);
  return 0;
}


/* Output:
size: 4
Exclusions: 2 3 4 5 
size: 4
Exclusions: 2 3 4 5 
size: 4
Exclusions: 2 3 4 5 
size: 3
Exclusions: 2 4 5 
 */
