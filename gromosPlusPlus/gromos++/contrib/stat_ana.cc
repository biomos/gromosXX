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
 * @file stat_ana.cc
 * calculate statistics from a set of values
 */

/**
 * @page programs Program Documentation
 *
 * @anchor stat_ana
 * @section stat_ana statistical analysis
 * @author @ref ns
 * @date 02.07.2009
 *
 * Program stat_ana can be used to calculate statistical quantities such as
 * averages, standard deviations, block averaged error estimate and distributions
 * from a given set of values.
 * The data file is specified by \@in and must contain an arbitrary number of lines
 * with a fixed number of columns of data per line.
 * If distributions are requested (\@dist) start and end values and the number of
 * bins have to be given for every column.
 *
 * Note that program @ref tcf can do the same kind of analysis.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@in</td><td>&lt;data file&gt; </td></tr>
 * <tr><td>[\@dist</td><td>&lt;start, end, bins for every column&gt;]</td></tr>
 * </table>
 *
 * Example:
 * @verbatim
   stat_ana
     @in         tser.out
     @dist       0 10 200  0 360 12
   @endverbatim

 * <hr>
 */

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/gmath/Stat.h"
#include "../src/gmath/Distribution.h"
#include "../src/gromos/Exception.h"

using namespace gmath;
using namespace args;
using namespace std;

void trim_comment(string & line) {
  size_t pos = line.find_first_of('#');
  if (pos != string::npos)
    line.erase(pos, line.length()-pos);
}

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "in" << "dist";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@in      <data file>\n";
  usage += "\t[@dist   <start end bins for every column>]\n";

  try {
    Arguments args(argc, argv, knowns, usage);
    ifstream in(args["in"].c_str());
    if (!in.is_open())
      throw gromos::Exception(argv[0], "Cannot open data file.");

    vector<Stat<double> > stat;
    bool do_dist = (args.count("dist")>=0);
    Arguments::const_iterator dist_iter = args.lower_bound("dist"),
            dist_end = args.upper_bound("dist");
    vector<Distribution> dist;
    unsigned int cols = 0;
    unsigned int line_num = 0;
    do {
      string line;
      getline(in, line);
      ++line_num;
      trim_comment(line);
      if (line.empty())
        continue;

      istringstream is(line);
      // in this for loop we try to detect the number of columns automatically
      for(double val; is >> val; ++cols) {
        stat.push_back(Stat<double>());
        stat[cols].addval(val);

        // create distributions?
        if (do_dist) {
          if (dist_iter == dist_end)
            throw gromos::Exception(argv[0], "Not enough @dist parameters given.");
          double start, end;
          unsigned int bins;
          if(!(istringstream(dist_iter++->second) >> start))
            throw gromos::Exception(argv[0], "In @dist, start is not numeric.");
          if(!(istringstream(dist_iter++->second) >> end))
            throw gromos::Exception(argv[0], "In @dist, end is not numeric.");
          if(!(istringstream(dist_iter++->second) >> bins))
            throw gromos::Exception(argv[0], "In @bins, start is not numeric.");
          dist.push_back(Distribution(start, end, bins));
          dist[cols].add(val);
        }
      }
      break;
    } while (!in.eof());

    if (cols == 0)
      throw gromos::Exception(argv[0], "No columns in data file found.");

    cerr << cols << endl;

    // read the rest of the file
    while(!in.eof()) {
      string line;
      getline(in, line);
      ++line_num;
      trim_comment(line);
      if (line.empty())
        continue;

      istringstream is(line);
      for(unsigned int i = 0; i < cols; ++i) {
        double val;
        is >> val;
        if (is.fail()) {
          ostringstream msg;
          msg << "Cannot read column " << i+1 << " in line " << line_num << ".";
          throw gromos::Exception(argv[0], msg.str());
        }

        stat[i].addval(val);
        if (do_dist)
          dist[i].add(val);

      } // for columns

    } // while lines

    in.close();

    // write averages etc.
    cout.precision(8);
    cout.setf(ios::left, ios::adjustfield);
    cout << setw(15) << "<x>:";
    cout.setf(ios::right, ios::adjustfield);
    for(unsigned int i = 0; i < cols; ++i) {
      cout << setw(15) << stat[i].ave();
    }
    cout << endl;
    cout.setf(ios::left, ios::adjustfield);
    cout << setw(15) << "ln<exp(x)>:";
    cout.setf(ios::right, ios::adjustfield);
    for(unsigned int i = 0; i < cols; ++i) {
      cout << setw(15) << stat[i].lnexpave();
    }

    cout << endl;
    cout.setf(ios::left, ios::adjustfield);
    cout << setw(15) << "rmsd:";
    cout.setf(ios::right, ios::adjustfield);
    for(unsigned int i = 0; i < cols; ++i) {
      cout << setw(15) << stat[i].rmsd();
    }
    cout << endl;
    cout.setf(ios::left, ios::adjustfield);
    cout << setw(15) << "ee:";
    cout.setf(ios::right, ios::adjustfield);
    for(unsigned int i = 0; i < cols; ++i) {
      cout << setw(15) << stat[i].ee();
    }
    cout << endl;
    cout.setf(ios::left, ios::adjustfield);
    cout << setw(15) << "n:";
    cout.setf(ios::right, ios::adjustfield);
    for(unsigned int i = 0; i < cols; ++i) {
      cout << setw(15) << stat[i].n();
    }
    cout << endl;
    cout.setf(ios::left, ios::adjustfield);
    cout << setw(15) << "min:";
    cout.setf(ios::right, ios::adjustfield);
    for(unsigned int i = 0; i < cols; ++i) {
      cout << setw(15) << stat[i].min();
    }
    cout << endl;
    cout.setf(ios::left, ios::adjustfield);
    cout << setw(15) << "max:";
    cout.setf(ios::right, ios::adjustfield);
    for(unsigned int i = 0; i < cols; ++i) {
      cout << setw(15) << stat[i].max();
    }
    cout << endl;
    
    if (do_dist) {
      for(unsigned int i = 0; i < cols; ++i) {
        ostringstream file_dist;
        file_dist << "dist_" << i+1 << ".dat";
        ofstream file(file_dist.str().c_str());
        dist[i].write(file);
        file.close();
        file_dist.clear();
        file_dist.str("");
        file_dist << "dist_norm_" << i+1 << ".dat";
        file.open(file_dist.str().c_str());
        dist[i].write_normalized(file);
        file.close();
      }
    }
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

