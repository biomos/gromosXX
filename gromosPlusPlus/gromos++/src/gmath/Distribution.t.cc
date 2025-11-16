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
#include "Distribution.h"

#include <iostream>

#include "../gromos/Exception.h"

using gmath::Distribution;

using namespace std;

int main() {
  try {

    Distribution di(0, 10, 5);
    cout << "distribution\n";
    di.write(cout);
    di.add(3.2);
    di.add(0.4);
    di.add(4.0);
    di.add(8.0);
    double test = 12;
    if (test != di.add(test)) cout << "value " << test
            << " out of range, not added\n";

    cout << "\ndistribution after adding elements\n";

    di.write(cout);
    cout << "average " << di.ave() << endl;
    cout << "rmsd " << di.rmsd() << endl;
    cout << "number of elements in section 4: " << di[3] << endl;
    cout << "middle value of section 4: " << di.value(3) << endl;

    return 0;
  }  catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}
