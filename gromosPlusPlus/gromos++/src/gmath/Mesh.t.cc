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
#include "Mesh.h"

#include <iostream>

#include "Vec.h"
#include "../gcore/Box.h"

using namespace std;
using namespace gcore;
using namespace gmath;

int main(int argc, char** argv) {
  Box b(Vec(10.0,  0.0,  0.0),
        Vec( 0.0, 10.0,  0.0),
        Vec( 2.0,  2.0, 10.0));

  Vec centre = b.K() + b.L() + b.M();
  centre *= 0.5;

  Mesh<double> m;
  m.setBox(b);
  m.setTitle("Nathans Mesh");
  m.resize(100, 60, 40);
  m = 0.0;

  for (int i = 0; i < m.size()[0]; ++i) {
    for (int j = 0; j < m.size()[1]; ++j) {
      for (int k = 0; k < m.size()[2]; ++k) {
        m(i, j, k) = (centre - m.pos(i, j, k)).abs();
      }
    }
  }

  m.write(cout);
  return 0;
}
