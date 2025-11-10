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
#include "Energy.h"

#include <cassert>
#include <iostream>
#include <cstdlib>
#include <string>

#include "../gio/InTopology.h"
#include "../gio/InG96.h"
#include "../gromos/Exception.h"
#include "PropertyContainer.h"
#include "../gcore/System.h"
#include "../gcore/GromosForceField.h"
#include "../bound/Boundary.h"
#include "../bound/RectBox.h"
#include "../gcore/MoleculeTopology.h"

using namespace gcore;
using namespace gio;
using namespace utils;
using namespace gmath;


using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << "Usage: " + string(argv[0]) + " <Topology> <coordinates>\n";
    exit(1);
  }
  try {
    InTopology it(argv[1]);
    System sys(it.system());
    GromosForceField gff(it.forceField());

    InG96 ic(argv[2]);
    ic.select("ALL");
    ic >> sys;
    ic.close();

    bound::Boundary *pbc;
    pbc = new bound::RectBox(&sys);

    /*
    string s="1:20";
    string t="a%1:1,2,3";
    cout << s << " " << t << endl;
    AtomSpecifier as(sys, s);
    PropertyContainer pc(sys, pbc);
    pc.addSpecifier(t);

    Energy en(sys, gff, *pbc);
    en.setAtoms(as);
    en.setProperties(pc);
    en.setCutOff(1.4);
    en.setRF(62.0, 0.0);
    en.calc();
    cout << en.vdw_m(0) << "\t" << en.el_m(0) << "\t"
         << en.vdw_s(0) << "\t" << en.el_s(0) << endl;
    cout << en.vdw(0) << "\t" << en.el(0) << endl;
    cout << en.vdw() << "\t" << en.el() << endl;
    cout << en.nb() << endl;
    cout << en.cov(0) << endl;
    cout << en.tot() << endl;
     */
    return 0;
  }  catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}

