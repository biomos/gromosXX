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

#include "SolventTopology.h"

#include <iostream>

#include "AtomTopology.h"
#include "Constraint.h"

using namespace gcore;
using namespace std;

int main() {
  SolventTopology mt;
  AtomTopology at;
  at.setIac(3);
  mt.addAtom(at);
  at.setIac(4);
  mt.addAtom(at);

  cout << "AtomTopology IAC: ";
  for (int i = 0; i < mt.numAtoms(); ++i)
    cout << mt.atom(i).iac() << ' ';
  cout << endl;

  Constraint b(1, 2);
  mt.addConstraint(b);
  b = Constraint(2, 3);
  b.setDist(0.1);

  mt.addConstraint(b);

  cout << "Constraints: ";
  ConstraintIterator bi(mt);
  for (; bi; ++bi)
    cout << '(' << bi()[0] << ',' << bi()[1] << ") " << bi().dist();
  cout << endl;

  mt.setSolvName("SPC");
  cout << "Solvname: ";
  cout << mt.solvName();
  cout << endl;

  return 0;
}
