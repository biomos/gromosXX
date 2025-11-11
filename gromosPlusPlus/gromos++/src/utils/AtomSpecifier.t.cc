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
#include "AtomSpecifier.h"

#include <cassert>
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>

#include "../gio/InTopology.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/Solvent.h"
#include "../gio/InG96.h"
#include "../gromos/Exception.h"

using namespace gcore;
using namespace gio;
using namespace utils;

using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 4) {
    cerr << "Usage: " + string(argv[0]) + " <Topology> <atomspecifier> <coordinates>\n";
    exit(1);
  }
  try {
    InTopology it(argv[1]);
    System sys(it.system());
    string s = argv[2];

    AtomSpecifier as(sys, s);

    AtomSpecifier bs(sys);
    string f = argv[3];
    gio::InG96 ic(f);
    cout << s << " consists of " << as.size() << " atoms:\n";
    for (int i = 0; i < as.size(); i++)
      cout << as.mol(i) << " : " << as.atom(i) << endl;
    ic.select("ALL");
    ic >> sys;
    ic.close();


    cout << "After reading coordinates  "
            << s << " consists of " << as.size() << " atoms:\n";

    as.sort();

    for (int i = 0; i < as.size(); i++)
      cout << as.mol(i) << " : " << as.atom(i) << endl;
    as.removeAtom(-1, 12);
    as.removeAtom(3);
    cout << "After removing something "
            << s << " consists of " << as.size() << " atoms:\n";
    sys.sol(0).setNumPos(0);

    cout << "After removing the solvent "
            << s << " consists of " << as.size() << " atoms:\n";


    cout << "\n\nto string: " << as.toString()[0] << "\n\n";

    cout << "positions:\n";
    for (int i = 0; i < as.size(); ++i) {
      cout << as.pos(i)[0] << endl;
    }

    /*    
    bs.addSpecifier("1:23-44");
    
    AtomSpecifier cs(sys);
    cs.addAtom(1,4);
    
    cs=as+bs;
    cout << endl;
    
    cout << endl <<  cs.size() << endl;
    for(int i=0; i< cs.size();i++)
      cout << cs.mol(i) << " : " << cs.atom(i) << endl;
    cs.removeAtom(1,1);
    
    cout << endl <<  cs.size() << endl;
    for(int i=0; i< cs.size();i++)
      cout << cs.mol(i) << " : " << cs.atom(i) << endl;

    cs.addAtom(1,6);
    cs.addType("CA");
    
    cout << endl <<  cs.size() << endl;
    for(int i=0; i< cs.size();i++)
      cout << cs.mol(i) << " : " << cs.atom(i) << endl;
     */



    return 0;
  }  catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}
