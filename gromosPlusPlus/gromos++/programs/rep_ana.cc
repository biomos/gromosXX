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
 * @page programs Program Documentation
 *
 * @anchor rep_ana
 * @section rep_ana analyze REMD output data
 * @author @ref mc cc
 * @date 3. 12. 2007
 *
 * Program rep_ana extracts the information stored in the replica
 * exchange molecular dynamics (REMD) output file replica.dat.
 * It produces seven multi column output files:
 * temperature.dat: run number vs. temperature of replica 
 * (column 2 corresponds to replica 1, column 3 to replica 2 etc.)
 * lambda.dat: run number vs. lambda value of replica
 * epot.dat: run number vs. potential energy of replica
 * probability.dat: run number vs. switching probability of replica
 * switches.dat: run number vs. switching data of replica 
 * (0 = no switch in this run, 1 = switch in this run)
 * prob_T.dat: switching probabilities per temperature
 * prob_l.dat: switching probabilities per lambda
 *
 * Furthermore it calculates an optimized temperature or lambda set
 * based on the fraction of replicas diffusing from the lowest
 * to the highest temperature (lambda value). This algorithm is 
 * based on:
 * Katzgraber, H. G.; Trebst, S.; Huse, D. A. & Troyer, M.
 * Feedback-optimized parallel tempering Monte Carlo
 * Journal of Statistical Mechanics-Theory And Experiment,
 * Iop Publishing Ltd, 2006, P03018"
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td>[\@repdata</td><td>&lt;REMD output file, replica.dat&gt;]</td></tr>
 * <tr><td>[\@traj</td><td>&lt;REMD master trajectory files&gt;]</td></tr>
 * <tr><td>[\@input</td><td>&lt;REMD input file&gt;]</td></tr>
 * </table>
 * Example:
 * @verbatim
 rep_ana
     @repdata       replica.dat
     @input         replica.imd
   @endverbatim

 * <hr>
 */

#include <cstdlib>
#include <iterator>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

#include "../src/args/Arguments.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace args;

#include "mk_script.h"

struct Replica_Data {
  int ID;
  int partner;
  unsigned int run;
  double Ti, Tj;
  double li, lj;
  double epot_i, epot_j;
  double p;
  int s;
};

int main(int argc, char *argv[]) {
  Argument_List knowns;
  knowns << "repdata";

  string usage = "# " + string(argv[0]);
  usage += "\n\t[@repdata      <REMD data file>]\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    vector<vector<double> > temp_series, lam_series, epot_series, switch_series, prob_series;
    unsigned int num_rep;
    vector<double> T;
    vector<double> lam;
    int num_T, num_l;

    ifstream rep(args["repdata"].c_str());

    // get number of temperatures
    string l, s;
    getline(rep, l);
    istringstream is(l);
    is >> s >> s >> s >> num_T;

    // get number of lambda values
    getline(rep, l);
    is.clear();
    is.str(l);
    is >> s >> s >> s >> s >> num_l;

    num_rep = num_T * num_l;

    // get list of temperatures
    getline(rep, l);
    is.clear();
    is.str(l);
    double d;
    is >> s;
    for (int i = 0; i < num_T; ++i) {
      is >> d;
      T.push_back(d);
    }

    // get list of lambda values
    getline(rep, l);
    is.clear();
    is.str(l);
    is >> s;
    for (int i = 0; i < num_l; ++i) {
      is >> d;
      lam.push_back(d);
    }

    getline(rep, l); // empty
    getline(rep, l); // header

    vector<double> entry;
    lam_series.resize(num_rep);
    temp_series.resize(num_rep);
    epot_series.resize(num_rep);
    switch_series.resize(num_rep);
    prob_series.resize(num_rep);
    for (unsigned int i = 0; i < num_rep; i++) {
      lam_series[i] = entry;
      temp_series[i] = entry;
    }
    unsigned int trials = 1;
    vector<double> old_temp(num_rep);
    vector<double> old_lam(num_rep);
    vector<double> new_temp(num_rep);
    vector<double> new_lam(num_rep);

    // initial set up
    for (unsigned int i = 0; i < num_rep; i++) {
      old_temp[i] = i + 1;
      old_lam[i] = i + 1;
    }
    
    int first_run = -1;

    // read in exchange data
    while (true) {
      getline(rep, l);
      if (rep.eof()) break;
      is.clear();
      is.str(l);

      Replica_Data r;
      is >> r.ID >> r.partner >> r.run >> r.li >> r.Ti >> r.epot_i
              >> r.lj >> r.Tj >> r.epot_j >> r.p >> r.s;

      // check for run
      if (r.run > trials) {
        ++trials;
        old_lam = new_lam;
        old_temp = new_temp;
      }
      // first run
      if (first_run == -1) {
        first_run = r.run;  // so this also works when user used equilibration steps and the first run != 1
      }

      if (int(r.run) == first_run) {
        // lambda
        if (r.s == 1) { // exchange successfull
          lam_series[r.ID - 1].push_back(r.li);
          lam_series[r.ID - 1].push_back(r.lj);
          new_lam[r.ID - 1] = r.partner;
        } else { // exchange failed
          lam_series[r.ID - 1].push_back(r.li);
          lam_series[r.ID - 1].push_back(r.li);
          new_lam[r.ID - 1] = r.ID;
        }
        // temperature
        if (r.s == 1) { // exchange successfull
          temp_series[r.ID - 1].push_back(r.Ti);
          temp_series[r.ID - 1].push_back(r.Tj);
          new_temp[r.ID - 1] = r.partner;
        } else { // exchange failed
          temp_series[r.ID - 1].push_back(r.Ti);
          temp_series[r.ID - 1].push_back(r.Ti);
          new_temp[r.ID - 1] = r.ID;
        }
        // epot, prob, switches
        epot_series[r.ID - 1].push_back(r.epot_i);
        epot_series[r.ID - 1].push_back(r.epot_j);
        prob_series[r.ID - 1].push_back(1);
        prob_series[r.ID - 1].push_back(r.p);
        switch_series[r.ID - 1].push_back(1);
        switch_series[r.ID - 1].push_back(r.s);
        
      } else { // other runs
        // lambda
        if (r.s == 1) { // exchange successfull
          for (unsigned int i = 0; i < num_rep; i++) {
            if (old_lam[i] == r.ID) {
              lam_series[i].push_back(r.lj);
              new_lam[i] = r.partner;
              break;
            }
          }
        } else { // exchange failed
          for (unsigned int i = 0; i < num_rep; i++) {
            if (old_lam[i] == r.ID) {
              lam_series[i].push_back(r.li);
              new_lam[i] = r.ID;
              break;
            }
          }
        }
        // temperature
        if (r.s == 1) { // exchange successfull
          for (unsigned int i = 0; i < num_rep; i++) {
            if (old_temp[i] == r.ID) {
              temp_series[i].push_back(r.Tj);
              new_temp[i] = r.partner;
              break;
            }
          }
        } else { // exchange failed
          for (unsigned int i = 0; i < num_rep; i++) {
            if (old_temp[i] == r.ID) {
              temp_series[i].push_back(r.Ti);
              new_temp[i] = r.ID;
              break;
            }
          }
        }
        // epot, prob, switch
        epot_series[r.ID - 1].push_back(r.epot_j);
        prob_series[r.ID - 1].push_back(r.p);
        switch_series[r.ID - 1].push_back(r.s);
      }
    } // while loop: read in

    cout << "# number of replicas: " << num_rep << endl;
    cout << "# read " << trials << " runs.\n";

    // write the files
    {
      ofstream temp("temperature.dat");
      temp.precision(4);
      temp.setf(ios::fixed, ios::floatfield);
      
      for (std::vector<vector<double> >::iterator it_trial = temp_series.begin(); it_trial != temp_series.end(); ++it_trial) {
        unsigned int i = std::distance( temp_series.begin(), it_trial ) + 1;
        temp << std::setw(10) << i;
        for (std::vector<double>::iterator it_num = it_trial->begin(); it_num != it_trial->end(); ++it_num) {
          temp << std::setw(12) << *it_num;
        }
        temp << "\n";
      }
      temp.close();
    }

    {
      ofstream lambda("lambda.dat");
      lambda.precision(4);
      lambda.setf(ios::fixed, ios::floatfield);

      for (std::vector<vector<double> >::iterator it_trial = lam_series.begin(); it_trial != lam_series.end(); ++it_trial) {
        unsigned int i = std::distance( lam_series.begin(), it_trial ) + 1;
        lambda << std::setw(10) << i;
        for (std::vector<double>::iterator it_num = it_trial->begin(); it_num != it_trial->end(); ++it_num) {
          lambda << std::setw(12) << *it_num;
        }
        lambda << "\n";
      }
      lambda.close();
    }

    {
      ofstream epot("epot.dat");
      epot.precision(4);
      epot.setf(ios::fixed, ios::floatfield);

      for (std::vector<vector<double> >::iterator it_trial = epot_series.begin(); it_trial != epot_series.end(); ++it_trial) {
        unsigned int i = std::distance( epot_series.begin(), it_trial ) + 1;
        epot << std::setw(10) << i;
        for (std::vector<double>::iterator it_num = it_trial->begin(); it_num != it_trial->end(); ++it_num) {
          epot << std::setw(18) << *it_num;
        }
        epot << "\n";
      }
      epot.close();
    }

    {
      ofstream prob("probability.dat");
      prob.precision(4);
      prob.setf(ios::fixed, ios::floatfield);

      for (std::vector<vector<double> >::iterator it_trial = prob_series.begin(); it_trial != prob_series.end(); ++it_trial) {
        unsigned int i = std::distance( prob_series.begin(), it_trial ) + 1;
        prob << std::setw(10) << i;
        for (std::vector<double>::iterator it_num = it_trial->begin(); it_num != it_trial->end(); ++it_num) {
          prob << std::setw(12) << *it_num;
        }
        prob << "\n";
      }
      prob.close();
    }

    {
      ofstream switches("switches.dat");
      switches.precision(4);
      switches.setf(ios::fixed, ios::floatfield);

      for (std::vector<vector<double> >::iterator it_trial = switch_series.begin(); it_trial != switch_series.end(); ++it_trial) {
        unsigned int i = std::distance( switch_series.begin(), it_trial ) + 1;
        switches << std::setw(10) << i;
        for (std::vector<double>::iterator it_num = it_trial->begin(); it_num != it_trial->end(); ++it_num) {
          switches << std::setw(12) << *it_num;
        }
        switches << "\n";
      }
      switches.close();
    }

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }

  return 0;
}
