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

#include "Solvent.h"

#include <cassert>
#include <iostream>

#include "SolventTopology.h"
#include "AtomTopology.h"
#include "Constraint.h"
#include "../gmath/Vec.h"

using namespace gcore;
using namespace std;
using namespace gmath;

using namespace std;

ostream & operator<<(ostream &o, Vec &v) {
  o << '(' << v[0] << ',' << v[1] << ',' << v[2] << ')';
  return o;
}

int main() {
  SolventTopology mt;
  AtomTopology at;
  at.setIac(3);
  mt.addAtom(at);
  at.setIac(4);
  mt.addAtom(at);

  Constraint b(1, 2);
  b.setDist(0.123);

  mt.addConstraint(b);
  b = Constraint(2, 3);
  mt.addConstraint(b);

  mt.setSolvName("SPC");

  Solvent solv(mt);

  cout << "Number of atoms in topo: " << solv.topology().numAtoms() << endl;
  cout << "Number of solvent coordinates: " << solv.numPos() << endl;

  solv.addPos(Vec(1, 2, 3));
  solv.addPos(Vec(4, 5, 6));
  solv.addPos(Vec(1, 2, 3));
  solv.addPos(Vec(4, 5, 6));
  solv.pos(2) = Vec(7, 8, 9);

  cout << "Number of solvent coordinates: " << solv.numPos() << endl;

  for (int i = 0; i < solv.numPos(); i++)
    cout << "Pos: " << i + 1 << " " << solv.pos(i) << endl;


  cout << "Constraints: ";
  ConstraintIterator bi(solv.topology());
  for (; bi; ++bi)
    cout << '(' << bi()[0] << ',' << bi()[1] << ") " << bi().dist();
  cout << endl;
  cout << "Number of solvent coordinates: " << solv.numPos() << endl;

  cout << "Number of solvent molecules: "
          << solv.numPos() / solv.topology().numAtoms() << endl;

  cout << "Setting number of Solv-Coords to 7" << endl;
  solv.setNumPos(7);
  cout << "Number of solvent Coords: ";
  cout << solv.numPos() << endl;

  return 0;
}
