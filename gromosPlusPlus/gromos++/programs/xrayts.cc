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

/**
 * @file xrayts.cc
 * x-ray time series
 */

/**
 * @page programs Program Documentation
 *
 * @anchor xrayts
 * @section xrayts x-ray time series
 * @author @ref ns
 * @date 14. 04. 09
 *
 * Extracts the crystallographic restraints information form a special trajectory.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@restraj</td><td>&lt;special trajectory file(s)&gt; </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
 xrayts
    @restraj  ex.trs
    @time     0 1
    @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <limits>
#include <string>

#include "../src/args/Arguments.h"
#include "../src/utils/groTime.h"
#include "../src/utils/RestrTraj.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace args;
using namespace utils;
using namespace gio;

int main(int argc, char** argv) {
  Argument_List knowns;

  knowns << "restraj" << "time";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@restraj    <restraint trajectory files>\n";
  usage += "\t[@time        <time dt>]\n";
  try {
    Arguments args(argc, argv, knowns, usage);
    Time time(args);
    RestrTraj trs;
    cout.precision(8);
    XrayRestrData min_inst, min_avg;
   
    min_inst.state().r_inst = min_avg.state().r_avg = numeric_limits<double>::max();
    // loop over all restraint trajectories
    for (Arguments::const_iterator
      iter = args.lower_bound("restraj"),
            to = args.upper_bound("restraj");
            iter != to; ++iter) {

      // open this restraint trajectory file
      trs.open((iter->second).c_str());
      // read in this restraint trajectory file
      while (!trs.eof()) {
        trs.read();
        XrayRestrData data;
        trs >> data >> time;
        cout << time
                << setw(15) << data.state().r_inst
                << setw(15) << data.state().r_free_inst
                << setw(15) << data.state().r_avg
                << setw(15) << data.state().r_free_avg << endl;
        if (min_inst.state().r_inst > data.state().r_inst) {
          min_inst = data;
        }
        if (min_avg.state().r_avg > data.state().r_avg) {
          min_avg = data;
        }
      }
    }
    cout << "# min. R inst: " << setw(15) << min_inst.state().r_inst
                << setw(15) << min_inst.state().r_free_inst
                << setw(15) << min_inst.state().r_avg
                << setw(15) << min_inst.state().r_free_avg << endl;
    cout << "# min. R avg: " << setw(15) << min_avg.state().r_inst
                << setw(15) << min_avg.state().r_free_inst
                << setw(15) << min_avg.state().r_avg
                << setw(15) << min_avg.state().r_free_avg << endl;

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

