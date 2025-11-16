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
 * @file eds_mult_all.cc
 * calculates the parameters s and EiR iteratively for multiple EDS states
 */

/**
 * @page programs Program Documentation
 *
 * @anchor eds_mult_all
 * @section eds_mult_all calculates the parameters s and EiR iteratively for multiple EDS states
 * @author @ref sr
 * @date 27. 01. 10
 *
 * Calculates iteratively the EDS parameters EiR and s from energy time series of
 * a number (numstat) of end states (vy) and the reference state R (vr).
 * The form of the used Hamiltonian (single s = 1, multiple s = 2, maximum
 * spanning tree = 3) has to be specified (form), as the number of s parameters
 * depends on the functional form. For form = 1, only a single s parameter is
 * calculated, for form = 2 N(N-1)/2 s parameters are calculated and for form = 3
 * (N-1) parameters, respectively. There are always N energy offset parameters EiR,
 * independent of the functional form. The same number of old parameters have
 * to be given (s and EiR).
 *
 * If a maximum spanning tree is used as functional form (form = 3), an initial tree must
 * be specified (tree) and if this tree shall be updated along with the
 * parameters (update_tree = 1) or not (update_tree = 0).
 *
 * The s parameters are calculated using
 * @f[
 * ln \sum_{j=1, j \neq i}^{M} \left[ \left( \left \langle e^{-\beta
 * (|\Delta V_{ji}| - \Delta E_{ji}^{R})} \right \rangle_{i}
 * \right)^{s} \right] = ln(M-1) - 1
 * @f]
 *
 * where M = N for form = 1, and M = 2 for form = 2,3, respectively.
 * The energy offset parameters are calculated using
 *
 * @f[
 * E_{i}^{R}(new) = - \beta^{-1} \cdot ln \left \langle \left( 1 +
 * \sum_{j=1, j \neq i}^{N} e^{-\beta (\Delta V_{ji} - \Delta E_{ji}^{R})}
 * \right)^{-1} \right \rangle_{R_{new}} + E_{i}^{R}(old)
 * @f]
 *
 * As the formulae are correlated, they are solved iteratively until both
 * parameters are converged.
 *
 * At the end of the program, the number of iterations are written out together
 * with the final parameters. For form = 3 and update_tree = 1, the new maximum spanning
 * tree is written to a separate file called tree.dat.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td>  \@temp</td><td>&lt;temperature for perturbation&gt; </td></tr>
 * <tr><td>  \@numstat</td><td>&lt;number of EDS states (= N)&gt; </td></tr>
 * <tr><td>  \@form</td><td>&lt;functional form (i.e. 1 s, N(N-1)/2 s or (N-1) s)&gt; </td></tr>
 * <tr><td>  \@vr</td><td>&lt;energy time series of state R&gt; </td></tr>
 * <tr><td>  \@vy</td><td>&lt;energy time series of states Y (N files)&gt; </td></tr>
 * <tr><td>  \@s</td><td>&lt;old s parameters (number according to functional form)&gt; </td></tr>
 * <tr><td>  \@EiR</td><td>&lt;old energy offset parameters (N values)&gt;</td></tr>
 * <tr><td>  [\@update_tree</td><td>&lt;new or updating of maximum spanning tree (required for form = 3)]&gt;</td></tr>
 * <tr><td>  [\@tree</td><td>&lt;maximum spanning tree (required for update_tree = 0 or 1)]&gt;</td></tr>
 * </table>
 *
 * Example 1 (with maximum spanning tree):
 * @verbatim
 eds_mult_all
    @temp       300
    @numstat      3
    @form         3
    @vr           eR.dat
    @vy           eA.dat eB.dat eC.dat
    @s            0.0324  0.0256
    @EiR          0.0  -5.4   4.3
    @update_tree  0 (no updating)
    @tree         tree.dat
    @endverbatim
 *
 * * Example 2 (with new maximum spanning tree):
 * @verbatim
 eds_mult_all
    @temp       300
    @numstat      3
    @form         3
    @vr           eR.dat
    @vy           eA.dat eB.dat eC.dat
    @s            0.0324  0.0256
    @EiR          0.0  -5.4   4.3
    @update_tree  2 (new tree)
    @endverbatim
 *
 * * Example 3:
 * @verbatim
 eds_mult_all
    @temp       300
    @numstat      3
    @form         2
    @vr           eR.dat
    @vy           eA.dat eB.dat eC.dat
    @s            0.0324  0.0256  0.0657
    @EiR          0.0  -5.4   4.3
    @endverbatim
 *
 * <hr>
 */

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <sstream>

#include "../src/args/Arguments.h"
#include "../src/gmath/Stat.h"
#include "../src/gmath/Physics.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace args;

/**
 * @struct pair_struct
 * pair struct to hold the data related to a EDS pair
 */
struct pair_struct {
  /**
   * To store the first pair member
   */
  unsigned int i;
  /**
   * To store the second pair member
   */
  unsigned int j;
  /**
   * To store the old s parameter
   */
  double s_old;
  /**
   * To store the new s parameter
   */
  double s_new;
  /**
   * To store the first energy offset
   */
  double EiR;
  /**
   * To store the second energy offset
   */
  double EjR;
  /**
   * To store the beta * dEiR
   */
  double term;
  /**
   * To store -ln<...>_i
   */
  double ln_i;
  /**
   * To store -ln<...>_j
   */
  double ln_j;
  /**
   * To store (V_B - V_A)
   */
  gmath::Stat<double> x;

  /**
   * constructor, sets everything to zero
   */
  pair_struct() : i(0), j(0), s_old(0.0), s_new(0.0), EiR(0.0), EjR(0.0),
  term(0.0), ln_i(0.0), ln_j(0.0) {
  }
};

gmath::Stat<double> read_data(string name, Arguments::const_iterator & iter);

void calculate_EiR(vector<pair_struct> & pair, gmath::Stat<double> & vr,
        vector<gmath::Stat<double> > & vy,
        vector<double> & eir_old, vector<double> & eir_new,
        double beta, int numpairs, int i,
        int numstat, double factor);

void calculate_s(pair_struct & pair);

void calculate_tree(vector<pair_struct> & pairs, vector<int> & treepairs,
                    int numstat);

int main(int argc, char** argv) {

  Argument_List knowns;

  knowns << "temp" << "numstat" << "form" << "vr" << "vy"
          << "s" << "EiR" << "tree" << "update_tree";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@temp          <temperature for perturbation>\n";
  usage += "\t@numstat        <number of EDS states (= N)>\n";
  usage += "\t@form           <functional form (i.e. 1 s, N(N-1)/2 s or (N-1) s)>\n";
  usage += "\t@vr             <energy time series of state R>\n";
  usage += "\t@vy             <energy time series of states Y (N files)>\n";
  usage += "\t@s              <old s parameters (number according to functional form)>\n";
  usage += "\t@EiR            <old energy offset parameters (N values)>\n";
  usage += "\t[@update_tree   <updating of maximum spanning tree>]\n";
  usage += "\t[@tree          <maximum spanning tree (required for @form = 3)>]\n";

  try {

    Arguments args(argc, argv, knowns, usage);

    // Get temperature
    double temp = args.getValue<double>("temp");
    const double beta = 1 / (gmath::physConst.get_boltzmann() * temp);

    // Get number of EDS states
    int numstat = args.getValue<int>("numstat");
    const int numpairs1 = numstat * (numstat - 1) / 2;
    const int numpairs3 = numstat - 1;

    // Get functional form
    unsigned int form = args.getValue<int>("form");
    vector<pair_struct> pairs;
    switch (form) {
      case 1:
        pairs.resize(numpairs1);
        break;
      case 2:
        pairs.resize(numpairs1);
        break;
      case 3:
        pairs.resize(numpairs3);
        break;
      default:
        throw gromos::Exception("eds_mult_all", "functional form must be 1,2 or 3");
        break;
    }

    // read the EiR
    args.check("EiR", 1);
    vector<double> eir_old(numstat);
    {
      int i = 0;
      for (Arguments::const_iterator iter = args.lower_bound("EiR"),
              to = args.upper_bound("EiR"); iter != to; ++iter, ++i) {
        std::istringstream is(iter->second.c_str());
        if (!(is >> eir_old[i])) {
          throw gromos::Exception("eds_mult_all", "EiR not numeric or too many values");
        }
      }
      if (i < numstat)
        throw gromos::Exception("eds_mult_all", "not enough EiR values (must be N)");
    }

    // read the s
    double smallest_s_old = 1.0;
    args.check("s", 1);
    {
      if (form > 1) {
        int i = 0;
        for (Arguments::const_iterator iter = args.lower_bound("s"),
                to = args.upper_bound("s"); iter != to; ++iter, ++i) {
          std::istringstream is(iter->second.c_str());
          if (!(is >> pairs[i].s_old)) {
            throw gromos::Exception("eds_mult_all", "s not numeric or too many values");
          }
        }
        if (form == 2 && i < numpairs1)
          throw gromos::Exception("eds_mult_all", "not enough s values (must be N(N-1)/2)");
        if (form == 3 && i < numpairs3)
          throw gromos::Exception("eds_mult_all", "not enough s values (must be (N-1))");
      } else {
        std::istringstream is(args["s"]);
        if (!(is >> smallest_s_old)) {
          throw gromos::Exception("eds_mult_all", "functional form not numeric");
        }
      }
    }

    // file to write final new tree to
    ofstream tree_out("tree_new.dat");

    // check update_tree if form = 3 and read tree
    unsigned int tree = 0;
    if (form == 3) {
      //update
      args.check("update_tree", 1);
      {
        std::istringstream is(args["update_tree"]);
        if (!(is >> tree)) {
          throw gromos::Exception("eds_mult_all", "updating tree switch not numeric");
        }
        if (tree != 0 && tree != 1)
          throw gromos::Exception("eds_mult_all", "update_tree must be 0 or 1");
      }
      // read tree
      if (tree == 0 || tree == 1) {
        args.check("tree", 1);
        ifstream file(args["tree"].c_str());
        if (!file.good()) {
          throw gromos::Exception("eds_mult_all", "Could not open file '" + args["tree"] + "'");
        }
        if (!file.is_open()) {
          throw gromos::Exception("eds_mult_all", "could not open file '" + args["tree"] + "'");
        }
        string line;
        stringstream linestream;
        int i = 0;
        int A, B;
        do {
          getline(file, line, '\n');
          if (!file.eof() && line[0] != '#') {
            linestream.clear();
            linestream.str(line);
            linestream >> A >> B;
            if (!linestream.good() && !linestream.eof()) {
              ostringstream os;
              os << "failed to read values from line\n" << line << "\ngot\n";
              throw gromos::Exception("eds_mult_all", os.str());
            }
            pairs[i].i = A - 1;
            pairs[i].j = B - 1;
            cout << "tree: pair = " << A << " and " << B << endl;
            ++i;
          } // if
        } while (!file.eof());
        file.close();
        if (i < numpairs3)
          throw gromos::Exception("eds_mult_all", "not enough pairs (must be (N-1))");
      }
    } else {
      // assign pairs
      unsigned int p = 0;
      for (int j = 0; j < numstat - 1; j++) {
        for (int k = j + 1; k < numstat; ++k, ++p) {
          pairs[p].i = j;
          pairs[p].j = k;
        }
      }
    }

    // read the reference state time series
    args.check("vr", 1);
    Arguments::const_iterator iter = args.lower_bound("vr");
    gmath::Stat<double> vr = read_data("vr", iter);

    // read the energy time series of all endstates
    vector<gmath::Stat<double> > vy(numstat);
    args.check("vy", 1);
    int i = 0;
    for (Arguments::const_iterator
      iter = args.lower_bound("vy"),
            to = args.upper_bound("vy");
            iter != to; ++iter, ++i) {
      gmath::Stat<double> vyi = read_data("vy", iter);
      vy[i] = vyi;
    }
    if (i < numstat)
      throw gromos::Exception("eds_mult_all", "not enough energy time series vy");

    //cout << "everything read in" << endl;

    // check whether time series have the same length
    if (vy[0].n() != vr.n())
      throw gromos::Exception("eds_mult_all", "Time series files differ in length!\n");

    //cout << "V_Y - V_R generated" << endl;
    cout.precision(6);
    cout << "# ############### New EDS Parameters: #################" << endl;

    unsigned int max_iter = 300;
    vector<double> eir_new(numstat);
    // convergence criterium
    double crit = 0.0000001;

    if (form == 1) {      
      double factor = 1.0 / (numstat - 1);

      // loop over pairs
      for (int j = 0; j < numpairs1; ++j) {
        // generate the energy difference time series from all endstates
        for (int k = 0; k < vy[pairs[j].j].n(); k++) {
          double diff = vy[pairs[j].j].data()[k] - vy[pairs[j].i].data()[k];
          pairs[j].x.addval(diff);
        }
        pairs[j].EiR = eir_old[pairs[j].i];
        pairs[j].EjR = eir_old[pairs[j].j];

        // calculate ln{<exp[-beta*|dVji|]>_i}
        gmath::Stat<double> weightA, weightB, statesumA, statesumB;
        // loop over data
        for (int k = 0; k < vr.n(); k++) {
          double temprefA = -beta*(vy[pairs[j].i].data()[k] - vr.data()[k]);
          double temprefB = -beta*(vy[pairs[j].j].data()[k] - vr.data()[k]);
          statesumA.addval(-beta*(abs(pairs[j].x.data()[k])) + temprefA);
          statesumB.addval(-beta*(abs(pairs[j].x.data()[k])) + temprefB);
          weightA.addval(temprefA);
          weightB.addval(temprefB);
        }
        pairs[j].ln_i = statesumA.lnexpave() - weightA.lnexpave();
        pairs[j].ln_j = statesumB.lnexpave() - weightB.lnexpave();
        
        pairs[j].s_old = smallest_s_old;
      }

      double diff_eir, diff_s;
      double smallest_s_new = 1.0;
      bool not_converged = true;
      unsigned int iterations = 0;
      do { // iterative procedure
        // calculate EiR: loop over states
        diff_eir = 0.0;
        for (int j = 0; j < numstat; ++j) {
          calculate_EiR(pairs, vr, vy, eir_old, eir_new,
                  beta, numpairs1, j, numstat, factor);
        }
        double E0R = eir_new[0];
        for (int j = 0; j < numstat; ++j) {
          // make EiR relative to first EiR
          eir_new[j] -= E0R;

          diff_eir += fabs(eir_new[j] - eir_old[j]);
          eir_old[j] = eir_new[j];
        }

        // calculate s: loop over pairs
        smallest_s_new = 1.0;
        for (int j = 0; j < numpairs1; ++j) {
          pairs[j].EiR = eir_old[pairs[j].i];
          pairs[j].EjR = eir_old[pairs[j].j];
          pairs[j].term = beta * (pairs[j].EjR - pairs[j].EiR);

          calculate_s(pairs[j]);

          if (pairs[j].s_new < smallest_s_new) {
            smallest_s_new = pairs[j].s_new;
          }
        }
        for (int j = 0; j < numpairs1; ++j) {
          pairs[j].s_old = smallest_s_new;
        }
        diff_s = fabs(smallest_s_new - smallest_s_old);
        smallest_s_old = smallest_s_new;

        if (diff_eir < crit && diff_s < crit) {
          not_converged = false;
        } else {
          iterations++;
          cout << "iteration: " << iterations << endl;
        }
      } while (not_converged && iterations < max_iter);
      cout << "# iterations: " << iterations << " of " << max_iter << "\n" << endl;

      // write final s parameters
      cout << "# new s parameter" << endl;
      cout << smallest_s_new << "\n" << endl;

    } else if (form == 2) {
      double factor = 1.0 / (numstat - 1);

      // loop over pairs
      for (int j = 0; j < numpairs1; ++j) {
        // generate the energy difference time series from all endstates
        for (int k = 0; k < vy[pairs[j].j].n(); k++) {
          double diff = vy[pairs[j].j].data()[k] - vy[pairs[j].i].data()[k];
          pairs[j].x.addval(diff);
        }
        pairs[j].EiR = eir_old[pairs[j].i];
        pairs[j].EjR = eir_old[pairs[j].j];

        // calculate ln{<exp[-beta*|dVji|]>_i}
        gmath::Stat<double> weightA, weightB, statesumA, statesumB;
        // loop over data
        for (int k = 0; k < vr.n(); k++) {
          double temprefA = -beta*(vy[pairs[j].i].data()[k] - vr.data()[k]);
          double temprefB = -beta*(vy[pairs[j].j].data()[k] - vr.data()[k]);
          statesumA.addval(-beta*(abs(pairs[j].x.data()[k])) + temprefA);
          statesumB.addval(-beta*(abs(pairs[j].x.data()[k])) + temprefB);
          weightA.addval(temprefA);
          weightB.addval(temprefB);
        }
        pairs[j].ln_i = statesumA.lnexpave() - weightA.lnexpave();
        pairs[j].ln_j = statesumB.lnexpave() - weightB.lnexpave();
      }

      double diff_eir, diff_s;
      bool not_converged = true;
      unsigned int iterations = 0;
      do { // iterative procedure
        // calculate EiR: loop over states
        diff_eir = 0.0;
        for (int j = 0; j < numstat; ++j) {
          calculate_EiR(pairs, vr, vy, eir_old, eir_new,
                  beta, numpairs1, j, numstat, factor);
        }
        double E0R = eir_new[0];
        for (int j = 0; j < numstat; ++j) {
          // make EiR relative to first EiR
          eir_new[j] -= E0R;

          diff_eir += fabs(eir_new[j] - eir_old[j]);
          eir_old[j] = eir_new[j];
        }

        // calculate s: loop over pairs
        diff_s = 0.0;
        for (int j = 0; j < numpairs1; ++j) {
          pairs[j].EiR = eir_old[pairs[j].i];
          pairs[j].EjR = eir_old[pairs[j].j];
          pairs[j].term = beta * (pairs[j].EjR - pairs[j].EiR);

          calculate_s(pairs[j]);

          diff_s += fabs(pairs[j].s_new - pairs[j].s_old);
          pairs[j].s_old = pairs[j].s_new;
        }

        if (diff_eir < crit && diff_s < crit) {
          not_converged = false;
        } else {
          iterations++;
          cout << "iteration: " << iterations << endl;
        }
      } while (not_converged && iterations < max_iter);
      cout << "# iterations: " << iterations << " of " << max_iter << "\n" << endl;

      // write final s parameters
      cout << "# new s parameters" << endl;
      for (int j = 0; j < numpairs1; ++j) {
        cout << pairs[j].s_new << " ";
      }
      cout << "\n" << endl;

    } else if (form == 3 && tree == 0) {
      double factor = double(numstat) / double(2 * (numstat - 1));

      // loop over pairs
      for (int j = 0; j < numpairs3; ++j) {
        // generate the energy difference time series from all endstates
        for (int k = 0; k < vy[pairs[j].j].n(); k++) {
          double diff = vy[pairs[j].j].data()[k] - vy[pairs[j].i].data()[k];
          pairs[j].x.addval(diff);
        }
        //pairs[j].x = x[j];
        pairs[j].EiR = eir_old[pairs[j].i];
        pairs[j].EjR = eir_old[pairs[j].j];

        // calculate ln{<exp[-beta*|dVji|]>_i}
        gmath::Stat<double> weightA, weightB, statesumA, statesumB;
        // loop over data
        for (int k = 0; k < vr.n(); k++) {
          double temprefA = -beta*(vy[pairs[j].i].data()[k] - vr.data()[k]);
          double temprefB = -beta*(vy[pairs[j].j].data()[k] - vr.data()[k]);
          statesumA.addval(-beta*(abs(pairs[j].x.data()[k])) + temprefA);
          statesumB.addval(-beta*(abs(pairs[j].x.data()[k])) + temprefB);
          weightA.addval(temprefA);
          weightB.addval(temprefB);
        }
        pairs[j].ln_i = statesumA.lnexpave() - weightA.lnexpave();
        pairs[j].ln_j = statesumB.lnexpave() - weightB.lnexpave();
      }

      double diff_eir, diff_s;
      bool not_converged = true;
      unsigned int iterations = 0;
      do { // iterative procedure
        // calculate EiR: loop over states
        diff_eir = 0.0;
        for (int j = 0; j < numstat; ++j) {
          calculate_EiR(pairs, vr, vy, eir_old, eir_new,
                  beta, numpairs3, j, numstat, factor);
        }
        double E0R = eir_new[0];
        for (int j = 0; j < numstat; ++j) {
          // make EiR relative to first EiR
          eir_new[j] -= E0R;
          diff_eir += fabs(eir_new[j] - eir_old[j]);
          eir_old[j] = eir_new[j];
          //cout << "E" << j << "R = " << eir_old[j] << endl;
        }

        // calculate s: loop over pairs
        diff_s = 0.0;
        for (int j = 0; j < numpairs3; ++j) {
          pairs[j].EiR = eir_old[pairs[j].i];
          pairs[j].EjR = eir_old[pairs[j].j];
          pairs[j].term = beta * (pairs[j].EjR - pairs[j].EiR);

          calculate_s(pairs[j]);

          diff_s += fabs(pairs[j].s_new - pairs[j].s_old);
          pairs[j].s_old = pairs[j].s_new;
          //cout << "s" << pairs[j].i << pairs[j].j << " = " << pairs[j].s_old << endl;
        }

        if (diff_eir < crit && diff_s < crit) {
          not_converged = false;
        } else {
          iterations++;
          cout << "iteration: " << iterations << endl;
        }
      } while (not_converged && iterations < max_iter);
      cout << "# iterations: " << iterations << " of " << max_iter << "\n" << endl;

      // write final s parameters
      cout << "# new s parameters" << endl;
      for (int j = 0; j < numpairs3; ++j) {
        cout << pairs[j].s_new << " ";
        tree_out << pairs[j].i + 1 << " " << pairs[j].j+ 1 << " " << pairs[j].s_new << endl;
      }
      cout << endl;
      tree_out.close();

    } else if (form == 3 && tree == 1) {
      double factor = double(numstat) / double(2 * (numstat - 1));

      // new pair_struct vector
      vector<pair_struct> allpairs(numpairs1);
      vector<int> treepairs(numpairs3);
      // assign pairs and s
      unsigned int p = 0, t = 0;
      for (unsigned int j = 0; int(j) < numstat - 1; j++) {
        for (unsigned int k = j + 1; int(k) < numstat; ++k, ++p) {
          allpairs[p].i = j;
          allpairs[p].j = k;
          allpairs[p].s_old = 1.0;
          for (unsigned int l = 0; int(l) < numpairs3; l++) {
            if (j == pairs[l].i && k == pairs[l].j) {
              treepairs[t] = p;
              allpairs[p].s_old = pairs[l].s_old;
              ++t;
            }
          }
        }
      }

      // loop over pairs
      for (int j = 0; j < numpairs1; ++j) {
        // generate the energy difference time series from all endstates
        for (int k = 0; k < vr.n(); k++) {
          double diff = vy[allpairs[j].j].data()[k] - vy[allpairs[j].i].data()[k];
          allpairs[j].x.addval(diff);
        }
        allpairs[j].EiR = eir_old[allpairs[j].i];
        allpairs[j].EjR = eir_old[allpairs[j].j];

        // calculate ln{<exp[-beta*|dVji|]>_i}
        gmath::Stat<double> weightA, weightB, statesumA, statesumB;
        // loop over data
        for (int k = 0; k < vr.n(); k++) {
          double temprefA = -beta*(vy[allpairs[j].i].data()[k] - vr.data()[k]);
          double temprefB = -beta*(vy[allpairs[j].j].data()[k] - vr.data()[k]);
          statesumA.addval(-beta*(abs(allpairs[j].x.data()[k])) + temprefA);
          statesumB.addval(-beta*(abs(allpairs[j].x.data()[k])) + temprefB);
          weightA.addval(temprefA);
          weightB.addval(temprefB);
        }
        allpairs[j].ln_i = statesumA.lnexpave() - weightA.lnexpave();
        allpairs[j].ln_j = statesumB.lnexpave() - weightB.lnexpave();
      }

      double diff_eir, diff_s;
      bool not_converged = true;
      unsigned int iterations = 0;
      do { // iterative procedure
        // update tree ?
        if (iterations == 1 || (iterations%(max_iter/10) == 0 && iterations != 0)) {
          cout << "update tree at iteration: " << iterations << endl;

          // now, get the new tree
          calculate_tree(allpairs, treepairs, numstat);

          // write the new tree back to the pairs vector
          for (int j = 0; j < numpairs3; ++j) {
            pairs[j] = allpairs[treepairs[j]];
          }
        }
        diff_eir = 0.0;
        for (int j = 0; j < numstat; ++j) {
          calculate_EiR(pairs, vr, vy, eir_old, eir_new,
                  beta, numpairs3, j, numstat, factor);
        }
        double E0R = eir_new[0];
        for (int j = 0; j < numstat; ++j) {
          // make EiR relative to first EiR
          eir_new[j] -= E0R;
          diff_eir += fabs(eir_new[j] - eir_old[j]);
          eir_old[j] = eir_new[j];
          //cout << "E" << j << "R = " << eir_old[j] << endl;
        }

        // calculate s: loop over pairs
        diff_s = 0.0;
        for (int j = 0; j < numpairs1; ++j) {
          allpairs[j].EiR = eir_old[allpairs[j].i];
          allpairs[j].EjR = eir_old[allpairs[j].j];
          allpairs[j].term = beta * (allpairs[j].EjR - allpairs[j].EiR);

          calculate_s(allpairs[j]);

          diff_s += fabs(allpairs[j].s_new - allpairs[j].s_old);
          allpairs[j].s_old = allpairs[j].s_new;
          //cout << "s" << allpairs[j].i << allpairs[j].j << " = " << allpairs[j].s_old << endl;
        }
        for (int j = 0; j < numpairs3; ++j) {
          pairs[j].s_old = allpairs[treepairs[j]].s_old;
          pairs[j].EiR = allpairs[treepairs[j]].EiR;
          pairs[j].EjR = allpairs[treepairs[j]].EjR;
        }

        if (diff_eir < crit && diff_s < crit) {
          not_converged = false;
        } else {
          iterations++;
          cout << "iteration: " << iterations << endl;
        }
      } while (not_converged && iterations < max_iter);
      cout << "# iterations: " << iterations << " of " << max_iter << "\n" << endl;

      // write final tree
      cout << "# new tree pairs" << endl;
      for (int j = 0; j < numpairs3; ++j) {
        cout << pairs[j].i + 1 << "  " << pairs[j].j + 1 << endl;
        tree_out << pairs[j].i + 1 << "  " << pairs[j].j + 1 << " " << pairs[j].s_old << endl;
      }
      tree_out.close();

      // write final s parameters
      cout << "# new s parameters" << endl;
      for (int j = 0; j < numpairs3; ++j) {
        cout << pairs[j].s_old << " ";
      }
      cout << endl;
    }


    // write final energy offsets
    cout << "# new energy offsets: " << endl;
    for (int j = 0; j < numstat; ++j) {
      cout << eir_old[j] << " ";
    }
    cout << "\n" << endl;

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }

  return 0;
}

gmath::Stat<double> read_data(string name, Arguments::const_iterator & iter) {

  gmath::Stat<double> data;

  // open the time series file for quantity x
  ifstream x;
  double time, q;

  x.open((iter->second).c_str());
  if (!x)
    throw gromos::Exception("eds_mult_all", "Could not open time series file for " + name + ".\n");

  // read
  string sdum;
  while (true) {
    std::getline(x, sdum);
    if (x.eof()) {
      break;
    }
    std::string::size_type it = sdum.find('#');
    if (it != std::string::npos)
      sdum = sdum.substr(0, it);
    if (sdum.find_first_not_of(" \t") == std::string::npos)
      continue;
    std::istringstream is(sdum);
    if ((!(is >> time >> q)) && !x.eof())
      throw gromos::Exception("eds_mult_all", "Error when reading from " + name + " time series file.\n");
    data.addval(q);

  }
  return data;
}

void calculate_EiR(vector<pair_struct> & pairs, gmath::Stat<double> & vr,
        vector<gmath::Stat<double> > & vy,
        vector<double> & eir_old, vector<double> & eir_new,
        double beta, int numpairs, int p_i,
        int numstat, double factor) {

  // calculate 1.0 / (sum_{j=1 and j != i}{exp(-beta*(dV_ji - dE_jiR))} + 1 )
  gmath::Stat<double> statesum, weight;
  for (int k = 0; k < vr.n(); ++k) {
    double sum_ln = 1.0;
    for (int j = 0; j < numstat; ++j) {
      if (j != p_i) {
          double deir = eir_old[j] - eir_old[p_i];
          double diff_ln = -beta * ((vy[j].data()[k]-vy[p_i].data()[k]) - deir);
          sum_ln = std::max(sum_ln, diff_ln)
              + log(1 + exp(std::min(sum_ln, diff_ln) - std::max(sum_ln, diff_ln)));
      }
    }
    // and V_Rnew
    double beta_s = -beta * pairs[0].s_old;
    double partA = beta_s * (vy[pairs[0].i].data()[k] - pairs[0].EiR);
    double partB = beta_s * (vy[pairs[0].j].data()[k] - pairs[0].EjR);
    double vrnew = (std::max(partA, partB)
                    + log(1 + exp(std::min(partA, partB) - std::max(partA, partB))))
                    / pairs[0].s_old;
    for (int i = 1; i < numpairs; ++i) {
      beta_s = -beta * pairs[i].s_old;
      partA = beta_s * (vy[pairs[i].i].data()[k] - pairs[i].EiR);
      partB = beta_s * (vy[pairs[i].j].data()[k] - pairs[i].EjR);
      double elem = (std::max(partA, partB)
                    + log(1 + exp(std::min(partA, partB) - std::max(partA, partB))))
                    / pairs[i].s_old;
      vrnew = std::max(vrnew, elem)
              + log(1 + exp(std::min(vrnew, elem) - std::max(vrnew, elem)));

    }
    vrnew = - (1.0 / beta) * (vrnew + log(factor));
    double diff = -beta * (vrnew - vr.data()[k]);
    statesum.addval(diff - sum_ln);
    weight.addval(diff);
  }

  // reweight
  double lnX_Y = statesum.lnexpave() - weight.lnexpave();

  // now calculate the new EiR
  eir_new[p_i] = (-(lnX_Y / beta) + eir_old[p_i]);

}

void calculate_s(pair_struct & pair) {
  long double sA = 1.0 / abs(-pair.ln_i - pair.term);
  long double sB = 1.0 / abs(-pair.ln_j + pair.term);
  pair.s_new = min(sA, sB);

  if (pair.s_new > 1)
    pair.s_new = 1;

}

void calculate_tree(vector<pair_struct> & pairs, vector<int> & treepairs,
                    int numstat) {
  vector<int> indicator(numstat, 1);
  indicator[0] = 0;
  unsigned int max_iter = 500;
  unsigned int iter = 0;
  bool not_converged = true;
  unsigned int p = 0;
  //for (unsigned int i = 0; i < max_iter; ++i) {
  do {
    int i_max = -5, j_max = -5;
    double s_max = -100.0;
    unsigned int mem_max = 0;
    for (int i = 0; i < numstat; ++i) {
      for (int j = 0; j < numstat; ++j) {
        if (j != i) {
          int mem = (i*(2*numstat - i - 1)/2) + j - i - 1;
          if (j < i)
            mem = (j*(2*numstat - j - 1)/2) + i - j - 1;
          if (s_max < pairs[mem].s_old && indicator[i] != indicator[j]) {
            s_max = pairs[mem].s_old;
            i_max = i;
            j_max = j;
            mem_max = mem;
          }
        }
      } // j loop
    } // i loop
    if (j_max > 0) {
      indicator[j_max] = 0;
      indicator[i_max] = 0;
      cout << "i_max = " << i_max << " , j_max = " << j_max << ", s_max = " << s_max << endl;
      assert(p < treepairs.size());
      treepairs[p] = mem_max;
      ++p;
      ++iter;
    } else {
      not_converged = false;
    }
  } while (not_converged && iter < max_iter); // iteration
  if (iter == max_iter)
    cout << "Warning: tree is not converged!" << endl;

}
