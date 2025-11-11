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
 * @file eds_update_2.cc
 * calculates the parameters s and EiR iteratively according to the procedures
 * described in J. Chem. Phys. 135 (2011) 024105 and J. Comput. Chem. 33 (2012) 640,
 * respectively.
 */

/**
 * @page programs Program Documentation
 *
 * @anchor eds_update_2
 * @section eds_update_2 calculates the parameters s and EiR iteratively
 * @author @ref nh
 * @date 11. 07. 12
 *
 * Calculates the EDS paramters EiR and s from energy time series of 
 * two endstates (vy) and the reference state R (vr). Two update 
 * schemes are implemented: Scheme 1 calculates the new parameters 
 * according to the procedure described in J. Chem. Phys. 135 (2011) 024105.
 * In that case the parameter eunder corresponds to the energy threshold.
 * Scheme 2 calculates new parameters according to the procedure 
 * described in J. Comput. Chem. 33 (2012) 640. In that case the 
 * parameter eunder corresponds to the energy separating sampling from 
 * state A from sampling of state B while the parameters etrans 
 * specifies the width of the transition region.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td>  \@temp</td><td>&lt;temperature for perturbation&gt; </td></tr>
 * <tr><td>  \@vr</td><td>&lt;energy time series of state R&gt; </td></tr>
 * <tr><td>  \@vy</td><td>&lt;energy time series of states Y (2 files)&gt; </td></tr>
 * <tr><td>  \@s</td><td>&lt;current s parameter&gt; </td></tr>
 * <tr><td>  \@s_old</td><td>&lt;old s parameter&gt; </td></tr>
 * <tr><td>  \@EiR</td><td>&lt;old energy offset parameters (2 values)&gt;</td></tr>
 * <tr><td>  \@update</td><td>&lt;which update scheme&gt;</td></tr>
 * <tr><td>  \@eunder</td><td>&lt;energy threshold if update=1; separation energy if update=2&gt;</td></tr>
 * <tr><td>  \@etrans</td><td>&lt;ignored if update=1; size of transition region if update=2&gt;</td></tr>
 * <tr><td>  [\@scale</td><td>&lt; scaling factor to modify default factors]&gt;</td></tr> 
 * </table>
 *
 * Example 1 (scheme according to J. Chem. Phys. 135 (2011) 024105):
 * @verbatim
 eds_update_2
    @temp       300
    @vr         eR.dat
    @vy         eA.dat eB.dat
    @s          0.0324
    @EiR        0.0  -5.4
    @update     1 
    @eunder     0.0
    @scale      1.0
    @endverbatim
 *
 * Example 2 (scheme according to J. Comput. Chem. 33 (2012) 640):
 * @verbatim
 eds_update_2
    @temp       300
    @vr         eR.dat
    @vy         eA.dat eB.dat 
    @s          0.0324  
    @EiR        0.0  -5.4 
    @update     2 
    @eunder     40.0
    @etrans     20.0
    @scale      1.0
    @endverbatim
 *
 * <hr>
 */

#include <algorithm>
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
   * To store the current s parameter
   */
  double s;
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
   * To store (V_B - V_A)
   */
  gmath::Stat<double> x;
  /**
   * To store the number of configurations belonging to the first state
   */
  unsigned int num_i;
  /**
   * To store the number of configurations belonging to the second state
   */
  unsigned int num_j;

  /**
   * constructor, sets everything to zero
   */
  pair_struct() : i(0), j(0), s_old(0.0), s(0.0), s_new(0.0), EiR(0.0), EjR(0.0),
  num_i(0), num_j(0){}
};

gmath::Stat<double> read_data(string name, Arguments::const_iterator & iter);

void calculate_EiR(vector<pair_struct> & pair, gmath::Stat<double> & vr,
        vector<gmath::Stat<double> > & vy,
        vector<double> & eir_old, vector<double> & eir_new,
        double beta, int numpairs, int i,
        int numstat, double factor);

void calculate_s(pair_struct & pair, vector<gmath::Stat<double> > & vy, 
        vector<double> & eir_old, int update, double eunder, double etrans, double scale);

int main(int argc, char** argv) {

  Argument_List knowns;

  knowns << "temp" << "vr" << "vy" << "s" << "s_old" << "EiR"
           << "update" << "eunder" << "etrans" << "scale";



  string usage = "# " + string(argv[0]);
  usage += "\n\t@temp          <temperature for perturbation>\n";
  usage += "\t@vr             <energy time series of state R>\n";
  usage += "\t@vy             <energy time series of states Y (2 files)>\n";
  usage += "\t@s              <current s parameter>\n";
  usage += "\t@s_old          <old s parameter>\n";
  usage += "\t@EiR            <energy offset parameters (2 values)>\n";
  usage += "\t@update         <which update scheme should be used>\n";
  usage += "\t@eunder         <energy threshold if update=1; separation energy if update=2>\n";
  usage += "\t@etrans         <ignored if update=1; size of transition region if update=2>\n";
  usage += "\t[@scale         <scaling factor to modify default factors>]\n";

  try {

    Arguments args(argc, argv, knowns, usage);

    // Get temperature
    double temp = args.getValue<double>("temp");
    const double beta = 1 / (gmath::physConst.get_boltzmann() * temp);

    // Get scale
    double scale = args.getValue<double>("scale", false, 1.0);

    // Get update scheme 
    int update = args.getValue<int>("update");
    if (!(update == 1 || update == 2)) {
      throw gromos::Exception("eds_update_2", "update must be 1 or 2");
    }

    // Get parameters for update scheme
    double eunder = args.getValue<double>("eunder");
    double etrans = 0.0;
    if (update == 2) {
      etrans = args.getValue<double>("etrans");
    }


    // Get number of EDS states
    int numstat = 2; // only two states allowed
    const int numpairs1 = numstat * (numstat - 1) / 2;

    vector<pair_struct> pairs;
    pairs.resize(numpairs1);

    // read the EiR
    args.check("EiR", 1);
    vector<double> eir_old(numstat);
    {
      int i = 0;
      for (Arguments::const_iterator iter = args.lower_bound("EiR"),
              to = args.upper_bound("EiR"); iter != to; ++iter, ++i) {
        std::istringstream is(iter->second.c_str());
        if (!(is >> eir_old[i])) {
          throw gromos::Exception("eds_update_2", "EiR not numeric or too many values");
        }
      }
      if (i < numstat)
        throw gromos::Exception("eds_update_2", "not enough EiR values (must be 2)");
    }

    // read the current s
    double smallest_s = 1.0;
    args.check("s", 1);
    {
      std::istringstream is(args["s"]);
      if (!(is >> smallest_s)) {
        throw gromos::Exception("eds_update_2", "current s parameter not numeric");
      }
    }

    // read the old s
    double smallest_s_old = 1.0;
    args.check("s_old", 1);
    {
      std::istringstream is(args["s_old"]);
      if (!(is >> smallest_s_old)) {
        throw gromos::Exception("eds_update_2", "old s parameter not numeric");
      }
    }     

    // assign pair
    unsigned int p = 0;
    for (int j = 0; j < numstat - 1; j++) {
      for (int k = j + 1; k < numstat; ++k, ++p) {
        pairs[p].i = j;
        pairs[p].j = k;
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
      throw gromos::Exception("eds_update_2", "not enough energy time series vy");

    //cout << "everything read in" << endl;

    // check whether time series have the same length
    if (vy[0].n() != vr.n())
      throw gromos::Exception("eds_update_2", "Time series files differ in length!\n");

    //cout << "V_Y - V_R generated" << endl;
    cout.precision(6);
    cout << "# ############### New EDS Parameters: #################" << endl;
    vector<double> eir_new(numstat);

    double factor = 1.0 / (numstat - 1);

    if (update == 1) {
    //double factor = 1.0 / (numstat - 1);

      //loop over states to calculate occupancies
      vector<unsigned int> occup(numstat);
      for (int j = 0; j < numstat; ++j) {
        occup[j] = 0;
      }
      for (int k = 0; k < vy[0].n(); ++k) {
        double enek = vy[0].data()[k];
        int state = 0;
        for (int j = 1; j < numstat; ++j) {
          if (vy[j].data()[k] < enek) {
            enek = vy[j].data()[k];
            state = j;
          }
        }
        occup[state]++;
      }
      for (int j = 0; j < numstat; ++j) {
        cout << "# state" << j << " = " << occup[j] << " ";
      }
      cout << endl;
    
      // loop over pairs
      for (int j = 0; j < numpairs1; ++j) {
        pairs[j].EiR = eir_old[pairs[j].i];
        pairs[j].EjR = eir_old[pairs[j].j];

        pairs[j].s_old = smallest_s_old;
        pairs[j].s = smallest_s;

        pairs[j].num_i = occup[pairs[j].i];
        pairs[j].num_j = occup[pairs[j].j];
      }

      // calculate s: loop over pairs
      double smallest_s_new = 1.0;
      for (int j = 0; j < numpairs1; ++j) {
        calculate_s(pairs[j], vy, eir_old, update, 
            eunder, etrans, scale);

        if (pairs[j].s_new < smallest_s_new) {
          smallest_s_new = pairs[j].s_new;
        }
      }

      for (int j = 0; j < numstat; ++j) {
        calculate_EiR(pairs, vr, vy, eir_old, eir_new,
                beta, numpairs1, j, numstat, factor);
      }
      double E0R = eir_new[0];
      for (int j = 0; j < numstat; ++j) {
        // make EiR relative to first EiR
        eir_new[j] -= E0R;
        eir_old[j] = eir_new[j];
      }

      // write final s parameters
      cout << "# new s parameter" << endl;
      cout << smallest_s_new << "\n" << endl;

    } else if (update == 2) {

      // loop over pairs
        for (int j = 0; j < numpairs1; ++j) {
          pairs[j].EiR = eir_old[pairs[j].i];
          pairs[j].EjR = eir_old[pairs[j].j];

          pairs[j].s_old = smallest_s_old;
          pairs[j].s = smallest_s;

         // pairs[j].num_i = occup[pairs[j].i];
         // pairs[j].num_j = occup[pairs[j].j];
        }

        double smallest_s_new = 1.0;
        for (int j = 0; j < numpairs1; ++j) {
          calculate_s(pairs[j], vy, eir_old, update, eunder, etrans, scale);

          if (pairs[j].s_new < smallest_s_new) {
            smallest_s_new = pairs[j].s_new;
          }
        }

        for (int j = 0; j < numstat; ++j) {
          calculate_EiR(pairs, vr, vy, eir_old, eir_new,
                  beta, numpairs1, j, numstat, factor);
        }
        double E0R = eir_new[0];
        for (int j = 0; j < numstat; ++j) {
          // make EiR relative to first EiR
          eir_new[j] -= E0R;
          eir_old[j] = eir_new[j];
        }

        // write final s parameters
        cout << "# new s parameter" << endl;
        cout << smallest_s_new << "\n" << endl;
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
    throw gromos::Exception("eds_update_2", "Could not open time series file for " + name + ".\n");

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
      throw gromos::Exception("eds_update_2", "Error when reading from " + name + " time series file.\n");
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
    double beta_s = -beta * pairs[0].s_new;
    double partA = beta_s * (vy[pairs[0].i].data()[k] - pairs[0].EiR);
    double partB = beta_s * (vy[pairs[0].j].data()[k] - pairs[0].EjR);
    double vrnew = (std::max(partA, partB)
                    + log(1 + exp(std::min(partA, partB) - std::max(partA, partB))))
                    / pairs[0].s_new;
    for (int i = 1; i < numpairs; ++i) {
      beta_s = -beta * pairs[i].s_new;
      partA = beta_s * (vy[pairs[i].i].data()[k] - pairs[i].EiR);
      partB = beta_s * (vy[pairs[i].j].data()[k] - pairs[i].EjR);
      double elem = (std::max(partA, partB)
                    + log(1 + exp(std::min(partA, partB) - std::max(partA, partB))))
                    / pairs[i].s_new;
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

void calculate_s(pair_struct & pair, vector<gmath::Stat<double> > & vy, vector<double> & eir_old, int update,
         double eunder, double etrans, double scale) {

// parameters for scheme 1:
  double downscale = 1.0 - ((scale*1.1-1.0)/2.0);
  double aveA = 0.0, aveB = 0.0;
  double ratio;

// parameters for scheme 2:
  int stateA = 0, stateB = 0, stateT = 0, numtrans = 0;
  int mintrans = 10;
  int visit = 0, visitold = 0;

  switch (update) {
    case 1:
      //double downscale = 1.0 - ((scale*1.1-1.0)/2.0);
      //double aveA = 0.0, aveB = 0.0;
      for (unsigned int i = 0; i < vy[pair.i].data().size(); i++) {
        aveA += vy[pair.i].data()[i];
        aveB += vy[pair.j].data()[i];
      }
      aveA /= vy[pair.i].data().size();
      aveB /= vy[pair.j].data().size();

      // up or down?
      //double ratio;
      if (aveA < eunder && aveB < eunder) { // always up!
        pair.s_new = pair.s * scale * 1.1;
        cout << "# s_old = " << pair.s_old << ", current s = " << pair.s << endl;
      } else {
      if (eir_old[pair.i] > eir_old[pair.j]) {
        if (pair.num_j > 0)
          ratio = double(pair.num_i) / double(pair.num_j);
        else
          ratio = pair.num_i;
        cout << "# s_old = " << pair.s_old << ", current s = " << pair.s << endl;
        if (ratio * 1.1 < 1 && ratio * 0.9 < 1 && pair.num_i != 0) { // up
          pair.s_new = pair.s * scale * 1.1;
        } else if ((ratio * 1.1 > 1 && ratio * 0.9 > 1) || pair.num_i == 0) { // down
          if (2*pair.s == pair.s_old) {
            pair.s_new = pair.s / 2;
          } else {
            pair.s_new = pair.s * downscale;
          }
        } else {
          pair.s_new = pair.s;
        }
      } else {
        if (pair.num_i > 0)
          ratio = double(pair.num_j) / double(pair.num_i);
        else
          ratio = pair.num_j;
        cout << "# s_old = " << pair.s_old << ", current s = " << pair.s << endl;
        if (ratio * 1.1 < 1 && ratio * 0.9 < 1 && pair.num_j != 0) { // up
          pair.s_new = pair.s * scale * 1.1;
        } else if ((ratio * 1.1 > 1 && ratio * 0.9 > 1) || pair.num_j == 0) { // down
          if (2*pair.s == pair.s_old) {
            pair.s_new = pair.s / 2;
          } else {
            pair.s_new = pair.s * downscale;
          }
        } else {
          pair.s_new = pair.s;
        }
      }
      }
      break;
//
    case 2:
      for (unsigned int i = 0; i < vy[pair.i].data().size(); i++) {
        if (vy[pair.j].data()[i] - vy[pair.i].data()[i] > eunder + etrans) {
          stateA++;
          visit = 0;
        } else if (vy[pair.j].data()[i] - vy[pair.i].data()[i] < eunder - etrans) {
          stateB++;
          visit = 1;
        } else {
          stateT++;
        }
        if (i > 0 && visitold != visit) {
          numtrans++;
        }
        visitold = visit;
      }  
 
      if (stateT > stateA && stateT > stateB) {
        pair.s_new = pair.s*1.1*scale;
      } else if (stateT >= stateA && numtrans > mintrans) { 
        pair.s_new = pair.s*1.05*scale;
      } else if (stateT >= stateB && numtrans > mintrans) { 
        pair.s_new = pair.s*1.05*scale;
      } else if ( (stateT < stateA || stateT < stateB ) && numtrans < mintrans ) {
        pair.s_new = pair.s*0.95*scale;
      } else if ( stateT*3.0 > stateA && stateT*3.0 > stateB && numtrans > mintrans ) {
        pair.s_new = pair.s*1.02*scale;
      } else {
        pair.s_new = pair.s;
      }
      cout << "# state A = " << stateA << ", state B = " << stateB << ", state T = " << stateT << endl;
      cout << "# number of A-B and B-A transitions = " << numtrans << endl;
      cout << "# s_old = " << pair.s_old << ", current s = " << pair.s << endl;
      break;
//
    default:
      throw gromos::Exception("eds_update_2","update must be 1 or 2");
      break;
  } // end of switch

}
