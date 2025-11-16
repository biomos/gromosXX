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
 * @file reweight.cc
 * reweight a time series
 */

/**
 * @page programs Program Documentation
 *
 * @anchor reweight
 * @section reweight reweight a time series
 * @author @ref cc
 * @date 04. 06. 09
 *
 * Reweights a time series of observed values of @f$X@f$ sampled during a simulation
 * at state @f$R@f$ (i.e. using the Hamiltonian @f$H_R=K_R(\vec{p})+V_R(\vec{r})@f$)
 * to another state @f$Y@f$ (neglecting kinetic contributions for simplicity):
 * @f[ \langle X \rangle_Y =
 *    \frac
 *        {\langle X \exp \left[-\beta \left (V_Y - V_R \right) \right] \rangle_R}
 *        {\langle   \exp \left[-\beta \left (V_Y - V_R \right) \right] \rangle_R}
 * = \langle X \exp \left[-\beta \left (V_Y - V_R -\Delta F_{YR} \right) \right] \rangle_R
 * @f]
 * with @f$\Delta F_{YR}=F_Y - F_R@f$.
 * The observed quantitiy @f$X@f$ can be a structural quantity (e.g. the time series
 * of an angle) or an energetic quantity (e.g. the time series of the ligand-protein
 * interaction energy). Note that the reweighting will only give useful results
 * if during the simulation at state @f$R@f$ all configurations that are important
 * to @f$Y@f$ are sampled.
 * The program reads three time series corresponding to the quantitiy @f$X@f$,
 * the energy of state @f$R@f$, and the energy of state @f$Y@f$. All time series
 * must have been calculated from the same ensemble @f$R@f$. The time series
 * files consist of a time column and a column containing the quantity (i.e. @f$X@f$,
 * @f$V_R@f$, or @f$V_Y@f$). The time series are obtained e.g. by @ref ene_ana or @ref tser.
 * If the bounds flag is given a normalized distribution of @f$X@f$ in the @f$Y@f$ ensemble will be
 * written out.
 * When calculating averages and distributions special care is taken
 * in order to avoid overflow (see Comput. Phys. Comm. 2003, 153, 397-406).
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td>  \@temp</td><td>&lt;temperature for perturbation&gt; </td></tr>
 * <tr><td>  \@x</td><td>&lt;time series of quantity X&gt; </td></tr>
 * <tr><td>  \@vr</td><td>&lt;energy time series of state R&gt; </td></tr>
 * <tr><td>  \@vy</td><td>&lt;energy time series of state Y&gt; </td></tr>
 * <tr><td> [\@bounds</td><td>&lt;lower bound&gt; &lt;upper bound&gt; &lt;grid points&gt;] </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
 reweight
    @temp      300
    @x         tser.dat
    @vr        eR.dat
    @vy        eY.dat
    @bounds    -300 300 100
    @endverbatim
 *
 * <hr>
 */

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <sstream>

#include "../src/args/Arguments.h"
#include "../src/gmath/WDistribution.h"
#include "../src/gmath/Physics.h"
#include "../src/gmath/Stat.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace args;

gmath::Stat<double> read_data(string name, Arguments &args);

int main(int argc, char** argv) {
  
  Argument_List knowns;

  knowns << "temp" << "x" << "vr" << "vy" << "bounds";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@temp     <temperature for perturbation>\n";
  usage +=   "\t[@x        <time series of quantity X>]\n";
  usage +=   "\t@vr       <energy time series of state R>\n";
  usage +=   "\t@vy       <energy time series of state Y>\n";
  usage +=   "\t[@bounds  <lower bound> <upper bound> <grid points>]\n";
 
  try {
    
    Arguments args(argc, argv, knowns, usage);

    // Get temperature as a double
    double temp = args.getValue<double>("temp", true);

    // read the bounds if given
    vector<double> bounds = args.getValues<double>("bounds", 3, false,
          Arguments::Default<double>() << 0.0 << 1.0 << 10.0);
    double dist_lower = bounds[0];
    double dist_upper = bounds[1];
    int dist_grid = int(bounds[2]);

    // read the time series
    bool has_x = false;
    gmath::Stat<double> x;
    if (args.count("x") > 0) {
      x = read_data("x", args);
      has_x = true;
    }
    gmath::Stat<double> vr = read_data("vr", args);
    gmath::Stat<double> vy =read_data("vy", args);

    // check whether all time series have the same length
    if ( x.n()!=vr.n() || x.n()!=vr.n()  )
      throw gromos::Exception("reweight", "Time series files differ in length!\n");

    // save -beta(V_Y - V_R)
    gmath::Stat<double> vyvr;
    gmath::Stat<double> Xexpvyvr;
    gmath::Stat<double> expvyvr;
    // create a distribution (with weights != 1)
    gmath::WDistribution xexpvyvr(dist_lower,dist_upper,dist_grid);
     
    /* loop over data that has been read in.
     * If speed turns out to be an issue this can be done
     * also on the fly when reading in the data
     */
    for (int i = 0; i < vr.n(); i++) {
      double diff = -(vy.data()[i] - vr.data()[i]) / (gmath::physConst.get_boltzmann() * temp);
      vyvr.addval(diff);
      expvyvr.addval(exp(diff));
      Xexpvyvr.addval(x.data()[i] * exp(diff));
      if (has_x)
        xexpvyvr.add(x.data()[i],diff);
    }
    
    if (has_x)
    

    cout.precision(10);
    // Calculate ln{|<X*exp[-beta(V_Y - V_R)]>_R|}
    int sign = 0;
    double lnXexpave = gmath::Stat<double>::lnXexpave(x, vyvr, sign);
    double lnexpave = vyvr.lnexpave();
    double dii = exp(2 * lnXexpave);
    double djj = exp(2 * lnexpave);
    double dji = exp(lnXexpave + lnexpave);
    
    // calculate statistical uncertainty
    double n = 1.0 / vr.n();
    double var_ii = gmath::Stat<double>::covariance(Xexpvyvr, Xexpvyvr);
    double si_ii = gmath::Stat<double>::stat_ineff(Xexpvyvr, Xexpvyvr);
    double d2i = var_ii * si_ii * n;
    double var_jj = gmath::Stat<double>::covariance(expvyvr, expvyvr);
    double si_jj = gmath::Stat<double>::stat_ineff(expvyvr, expvyvr);
    double d2j = var_jj * si_jj * n;
    double var_ji = gmath::Stat<double>::covariance(Xexpvyvr, expvyvr);
    double si_ji = gmath::Stat<double>::stat_ineff(Xexpvyvr, expvyvr);
    double d2ji = var_ji * si_ji * n;
    double error = sqrt( d2i/dii + d2j/djj - 2*d2ji/dji );
    
    cout << "# ln{|<X*exp[-beta(V_Y - V_R)]>_R|} = " << lnXexpave << endl;
    cout << "# sign = " << sign << endl;
    // Calculate ln{<exp[-beta(V_Y - V_R)]>_R}
    cout << "# ln{<exp[-beta(V_Y - V_R)]>_R} = " << lnexpave << endl;
    // <X>_Y
    cout << "# <X>_Y = " << exp(lnXexpave - lnexpave) * sign << setw(18) << error << endl;

    // Write out a distribution if the @bounds flag is given
    if (args.count("bounds") >= 0) {
      // write a normalized distribution
      xexpvyvr.write_normalized(cout);
      //xexpvyvr.write(cout);
    }
  
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
   
  return 0;
}

gmath::Stat<double> read_data(string name, Arguments & args) {

  gmath::Stat<double> data;

  // make sure we have the @name flag
  args.check(name.c_str(), 1);

  // open the time series file for quantity x
  ifstream x;
  double time, q;
  Arguments::const_iterator iter2 = args.lower_bound(name.c_str());

  x.open((iter2->second).c_str());
  if (!x)
    throw gromos::Exception("reweight", "Could not open time series file for " + name + ".\n");

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
      throw gromos::Exception("reweight", "Error when reading from " + name + " time series file.\n");
    data.addval(q);
    
  }
  return data;
}
