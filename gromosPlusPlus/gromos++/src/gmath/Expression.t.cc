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
#include "Expression.h"

#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>

#include "../gromos/Exception.h"

using gmath::Expression;

using namespace std;

int main(int argc, char *argv[]) {
  try {
    if (argc < 2) {
      cerr << "Usage: " + string(argv[0]) + " <expression>" << endl;
      exit(1);
    }
    ostringstream os;
    for (int i = 1; i < argc; i++)
      os << argv[i] << " ";

    string s = os.str();
    Expression e(s);
    e.writeExpression(cout);
    cout << e.value() << endl;
  }  catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}
