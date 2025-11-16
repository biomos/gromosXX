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
 * @file tcf.cc
 * calculate distributions and time correlation functions
 */

/**
 * @page programs Program Documentation
 *
 * @anchor tcf
 * @section tcf calculate distributions and time correlation functions
 * @author @ref co
 * @date 28. 7. 2006
 *
 * Program tcf performs simple statistical analyses, calculates distributions 
 * and time-correlation functions for any series of data points. As input it 
 * takes files with data listed in any number of columns but no further 
 * formatting, e.g. the output files of programs @ref tser or @ref ene_ana . 
 * Lines starting with the character "#" are ignored.
 *
 * For data in the specified columns, the program writes out the number of data
 * points, the average value, root-mean-square fluctuations, a statistical 
 * error estimate as well as the minimal and maximal value observed. The error 
 * estimate is calculated from block averages of different sizes, as described 
 * in Allen and Tildesley: "Computer Simulation of Liquids", 1987. In addition 
 * the program can calculate distributions and time-correlation functions. The 
 * program does not read time from the data-file, but the time interval between
 * data points can be specified by the user. Otherwise it is taken to be 1. 
 *
 * Distributions can be calculated for data in specified columns and can be 
 * normalized. 
 *
 * Time correlation functions of the general form 
 * @f[ C(t) = \left< f(A(\tau),B(\tau+t)) \right>_{\tau} @f]
 * can be calculated, where @f$A(\tau)@f$ and @f$B(\tau+t)@f$ represent the
 * data points at different time points and the user can specify any function
 * f(A,B). The program can calculate both auto-correlation functions (B=A) and 
 * cross correlation functions (B!=A) for time series of scalars or vectors. 
 * If A and B are represented by scalars and 
 * @f$f(A(\tau),B(\tau+t)) = A(\tau) * B(\tau+t)@f$, the program makes use of 
 * fast fourier transforms to calculate C(t). In other cases a direct summing 
 * algorithm is used, which may be considerably slower.
 *
 * In cases where one is interested in the correlation function of the 
 * fluctuations around the average, this average value can be subtracted from 
 * the data points before calculating the correlation function. A power 
 * spectrum can also be calculated. Because the correlation function is usually
 * very noisy at larger times, the noise level can be specified as the fraction
 * of the correlation function that should be considered. This part of the 
 * correlation function is then smoothened by multiplication by a cosine to 
 * make sure that it is zero at the end. It is then mirrored: all data points 
 * are repeated in reverse order at the end. From this the fourier transform is
 * taken, which is the resulting spectrum.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@files</td><td>&lt;data files&gt; </td></tr>
 * <tr><td> \@time</td><td>&lt;@ref utils::Time "time and dt"&gt; </td></tr>
 * <tr><td> [\@distribution</td><td>&lt;data columns to consider&gt;] </td></tr>
 * <tr><td> [\@bounds</td><td>&lt;lower bound&gt; &lt;upper bound&gt; &lt;grid points&gt;] </td></tr>
 * <tr><td> [\@normalize</td><td>(normalize the distributions)] </td></tr>
 * <tr><td> [\@tcf</td><td>&lt;data columns to consider&gt;] </td></tr>
 * <tr><td> [\@expression</td><td>&lt;expression for correlation function&gt;] </td></tr>
 * <tr><td> [\@spectrum</td><td>&lt;noise level&gt;] </td></tr>
 * <tr><td> [\@subtract_average</td><td>(take difference with respect to average value for tcf)] </td></tr>
 * </table>
 *
 * <b>See also:</b> @ref gmath::Correlation , @ref gmath::Distribution ,
 * @ref gmath::Expression , @ref gmath::Stat
 *
 * Example:
 * @verbatim
  tcf
    @files          tser.out
    @time           0 0.5
    @distribution   2
    @bounds         -180 180 20
    @normalize    
    @tcf            2
    @expression     "3 * cos ( A * B ) - 1" 
#   @spectrum
#   @subtract_average
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <set>

#include "../src/args/Arguments.h"
#include "../src/gmath/Stat.h"
#include "../src/gmath/Distribution.h"
#include "../src/gmath/Correlation.h"
#include "../src/gmath/Vec.h"
#include "../src/gromos/Exception.h"

using namespace args;
using namespace std;

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "files" << "distribution" << "normalize" << "bounds" << "tcf"
          << "expression" << "spectrum" << "subtract_average" << "time";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@files              <data files>\n";
  usage += "\t@time               <time> <time step>\n";
  usage += "\t[@distribution      <data columns to consider>]\n";
  usage += "\t[@bounds            <lower bound> <upper bound> <grid points>]\n";
  usage += "\t[@normalize         (normalize the distributions)]\n";
  usage += "\t[@tcf               <data columns to consider>]\n";
  usage += "\t[@expression        <expression for correlation function>]\n";
  usage += "\t[@spectrum          <noise level>]\n";
  usage += "\t[@subtract_average  (take difference with respect to average value for tcf)]\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // read the time
    vector<double> timearg = args.getValues<double>("time", 2, false,
            Arguments::Default<double>() << 0.0 << 1.0);
    double time = timearg[0];
    double dt = timearg[1];

    // determine which data to store
    set<int> data_sets;
    vector<int> data_index;

    // get all input for distributions
    bool distribution = false;
    vector<int> dist_index;
    for (Arguments::const_iterator iter = args.lower_bound("distribution"),
            to = args.upper_bound("distribution"); iter != to; ++iter) {
      int i = atoi(iter->second.c_str()) - 1;
      dist_index.push_back(i);
      data_sets.insert(i);
      distribution = true;
    }

    bool dist_normalize = false;
    if(args.count("normalize")>=0) dist_normalize=true;

    vector<double> bounds = args.getValues<double>("bounds", 3, false,
            Arguments::Default<double>() << 0.0 << 1.0 << 10.0);
    double dist_lower = bounds[0];
    double dist_upper = bounds[1];
    int dist_grid = int(bounds[2]);

    bool tcf = false;
    bool tcf_vector = false;
    bool tcf_expression = false;
    bool tcf_spectrum = false;
    bool tcf_subtract = false;
    bool tcf_auto = false;

    string tcf_expression_string;
    double tcf_noise = 1.0;
    vector<int> tcf_index;
    for (Arguments::const_iterator iter = args.lower_bound("tcf"),
            to = args.upper_bound("tcf"); iter != to; ++iter) {
      int i = atoi(iter->second.c_str()) - 1;
      tcf_index.push_back(i);
      data_sets.insert(i);
      tcf = true;
    }

    // handle the cases for autocorrelation functions
    // and vector correlation functions
    if (tcf) {
      if (tcf_index.size() == 1) {
        tcf_index.push_back(tcf_index[0]);
        tcf_auto = true;
      } else if (tcf_index.size() == 3) {
        tcf_index.push_back(tcf_index[0]);
        tcf_index.push_back(tcf_index[1]);
        tcf_index.push_back(tcf_index[2]);
        tcf_vector = true;
        tcf_auto = true;
      } else if (tcf_index.size() == 6) tcf_vector = true;
      else if (tcf_index.size() != 2)
        throw gromos::Exception("protcf",
              "Specify either 1, 2, 3, or 6 columns for tcf");
    }

    if (args.count("expression") > 1) {
      tcf_expression = true;
      for (Arguments::const_iterator iter = args.lower_bound("expression"),
              to = args.upper_bound("expression"); iter != to; ++iter) {
        tcf_expression_string += iter->second + " ";
      }
    }
    if (args.count("subtract_average") >= 0) tcf_subtract = true;
    if (args.count("spectrum") >= 0) {
      tcf_spectrum = true;
      if (args.count("spectrum") == 1) tcf_noise = atof(args["spectrum"].c_str());
    }

    // process the data sets that need to be stored.
    int data_max = 0;
    for (set<int>::const_iterator iter = data_sets.begin(),
            to = data_sets.end();
            iter != to; ++iter) {

      bool keep = true;
      if (tcf_vector) {
        keep = false;
        for (unsigned int i = 0; i < dist_index.size(); i++)
          if (*iter == dist_index[i]) keep = true;
      }
      if (keep)
        data_index.push_back(*iter);
      if (*iter > data_max) data_max = *iter;
    }
    data_max++;

    // calculate sizes for datastructures
    int data_num = data_index.size();

    // we will also need the inverted data_index
    vector<int> data_inv(data_max);
    for (int i = 0; i < data_num; i++)
      data_inv[data_index[i]] = i;

    // create data structures to read data into
    vector<gmath::Stat<double> > data;
    for (int i = 0; i < data_num; i++) {
      gmath::Stat<double> tmp;
      data.push_back(tmp);
    }

    vector<vector<gmath::Vec> > data_vec(2);
    vector<double> tmp_data(data_max);

    if (args.count("files") <= 0)
      throw gromos::Exception(argv[0], "There is no data file specified\n" +
            usage);

    // ALL INPUT GATHERED:
    // loop over the files and store data
    string line;
    stringstream linestream;
    for (Arguments::const_iterator
      iter = args.lower_bound("files"),
            to = args.upper_bound("files");
            iter != to; ++iter) {

      ifstream file(iter->second.c_str());
      if (!file.good()) {
        throw gromos::Exception("tcf", "Could not open file '" + iter->second + "'");
      }
      if (!file.is_open()) {
        throw gromos::Exception("tcf", "could not open file '" + iter->second + "'");
      }

      do {
        getline(file, line, '\n');
        if (!file.eof() && line[0] != '#') {
          linestream.clear();
          linestream.str(line);
          for (int i = 0; i < data_max; i++)
            linestream >> tmp_data[i];
          if (!linestream.good() && !linestream.eof()) {
            ostringstream os;
            os << "failed to read " << data_max << " values from line\n"
                    << line << "\ngot\n";
            for (int i = 0; i < data_max; i++)
              os << tmp_data[i] << "  ";
            throw gromos::Exception("protcf", os.str());
          }

          for (unsigned int j = 0; j < data_index.size(); j++)
            data[j].addval(tmp_data[data_index[j]]);
          if (tcf_vector) {
            gmath::Vec v1(tmp_data[tcf_index[0]],
                    tmp_data[tcf_index[1]],
                    tmp_data[tcf_index[2]]);
            data_vec[0].push_back(v1);
            gmath::Vec v2(tmp_data[tcf_index[3]],
                    tmp_data[tcf_index[4]],
                    tmp_data[tcf_index[5]]);
            data_vec[1].push_back(v2);
          }
        }
      } while (!file.eof());


      file.close();
    }

    // OK, we have all the data, now we need to do whatever we are expected
    // to do

    // CALCULATIONS AND OUTPUT
    cout << "TITLE\n"
            << "Statistical analysis of data file";
    if (args.count("files") > 1) cout << "s";
    cout << ":\n";
    for (Arguments::const_iterator
      iter = args.lower_bound("files"),
            to = args.upper_bound("files");
            iter != to; ++iter) {
      cout << iter->second << "\n";
    }
    cout << "END\n";

    // averages etc.
    if (data_index.size()) {
      cout << "STATISTICS\n"
              << "# column       N     average        rmsd  error est."
              << "     minimum     maximum\n";
      for (unsigned int i = 0; i < data_index.size(); i++) {
        int di = data_inv[data_index[i]];
        cout << setw(8) << data_index[i] + 1
                << ' ' << setw(7) << data[di].n()
                << ' ' << setw(11) << data[di].ave()
                << ' ' << setw(11) << data[di].rmsd()
                << ' ' << setw(11) << data[di].ee()
                << ' ' << setw(11) << data[di].min()
                << ' ' << setw(11) << data[di].max()
                << "\n";
      }
      cout << "END\n";
    }

    // create the distributions
    if (distribution) {
      cout << "DISTRIBUTION\n"
              << "# distributions calculated for column";
      if (dist_index.size() > 1) cout << "s";
      for (unsigned int i = 0; i < dist_index.size(); i++) {
        int di = data_inv[dist_index[i]];
        data[di].dist_init(dist_lower, dist_upper, dist_grid);
        cout << " " << dist_index[i] + 1;
      }
      cout << "\n# lower bound: " << dist_lower
              << "\n# upper bound: " << dist_upper
              << "\n# number of grid points: " << dist_grid
              << "\n# distribution is ";
      if (!dist_normalize) cout << "not ";
      cout << "normalized\n\n"
              << "#    value";
      for (unsigned int i = 0; i < dist_index.size(); i++)
        cout << setw(6) << dist_index[i] + 1 << ". column";
      cout << "\n";

      for (int i = 0; i < dist_grid; i++) {
        int dd = data_inv[dist_index[0]];
        cout << setw(10) << data[dd].distribution().value(i);
        for (unsigned int j = 0; j < dist_index.size(); j++) {
          int dj = data_inv[dist_index[j]];
          if (dist_normalize)
            cout << setw(14) << double((data[dj].distribution())[i]) /
            (data[dj].distribution().nVal()*
                  (data[dj].distribution().value(1) -
                  data[dj].distribution().value(0)));
          else
            cout << setw(14) << (data[dj].distribution())[i];
        }
        cout << "\n";
      }
      cout << "END\n";
    }
    if (tcf) {
      cout << "TIME CORRELATION FUNCTION\n"
              << "# calculating ";
      if (tcf_auto) cout << "auto-";
      cout << "correlation function for ";
      if (tcf_vector) {
        cout << "vector";
        if (!tcf_auto) cout << "s";
        cout << " defined by ";
      }
      cout << "column";
      if (!tcf_auto || tcf_vector) cout << "s";
      int limit = tcf_index.size();
      if (tcf_auto) limit /= 2;
      for (int i = 0; i < limit; i++) {
        if (i == int(tcf_index.size()) / 2) cout << " (A) and";
        cout << " " << tcf_index[i] + 1;
      }
      if (!tcf_auto) cout << " (B)";
      cout << "\n\n";
      if (tcf_subtract)
        cout << "# average values are subtracted from time series\n\n";
      cout << "# correlation function calculated as C(t) = <";
      if (tcf_expression) {
        cout << " f( A(T), ";
        if (tcf_auto) cout << "A";
        else cout << "B";
        cout << "(T+t) ) ";
      } else {
        cout << " A(T) * ";
        if (tcf_auto) cout << "A";
        else cout << "B";
        cout << "(T+t) ";
      }
      cout << ">_T\n";
      if (tcf_expression) {
        cout << "# with f( A(T), ";
        if (tcf_auto) cout << "A";
        else cout << "B";
        cout << "(T+t) ) = " << tcf_expression_string << "\n";
      }
      cout << "# using ";
      if (tcf_expression || tcf_vector)
        cout << "a double loop algorithm\n";
      else
        cout << "fast fourier transforms\n";

      gmath::Correlation *corr;
      // several cases
      if (tcf_vector) {
        corr = new gmath::Correlation(data_vec[0], data_vec[1]);
        corr->calc_direct();
      } else {
        int d1 = data_inv[tcf_index[0]];
        int d2 = data_inv[tcf_index[1]];
        if (tcf_subtract) {
          data[d1].subtract_average();
          data[d2].subtract_average();
        }
        corr = new gmath::Correlation(data[d1], data[d2]);
        if (tcf_expression) {
          corr->calc_expression(tcf_expression_string);
        }
        else {
          corr->calc_fft();
        }

      }
      cout << "\n#        t          C(t)\n";

      for (unsigned int i = 0; i < corr->size(); i++, time += dt) {
        cout << setw(10) << time
                << setw(14) << (*corr)[i] << "\n";
      }
      cout << "END\n";

      if (tcf_spectrum) {
        cout << "SPECTRUM\n"
                << "# calculated from above correlation function\n"
                << "# " << tcf_noise * 100.0 << " % of C(t) used in spectrum "
                << "calculation\n\n"
                << "#  frequency   intensity\n";

        vector<double> freq, spec;
        corr->spectrum(freq, spec, dt, tcf_noise);
        for (unsigned int i = 0; i < freq.size(); i++) {
          cout << setw(12) << freq[i] << setw(14) << spec[i] << "\n";
        }
        cout << "END\n";
      }

    }


  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}




