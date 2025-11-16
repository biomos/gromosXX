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

// args_Arguments.t.cc

#include "Arguments.h"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>

#include "../gromos/Exception.h"

using namespace std;
using namespace args;

int debug_level = 0;

namespace std {

  std::ostream & operator<<(std::ostream &os, const Arguments &args) {
    for (Arguments::const_iterator iter = args.begin();
            iter != args.end(); ++iter) {
      os << iter->first << ' ' << iter->second << endl;
    }
    return os;
  }
}

int main(int argc, char **argv) {
  Argument_List knowns;
  knowns << "bla" << "gug";
  string usage = argv[0];
  usage += " @bla <testargs> @gug <testargs>";

  try {
    Arguments args(argc, argv, knowns, usage);
    cout << args.count("gug") << endl;

    cout << args;


  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }

  return 0;
}
