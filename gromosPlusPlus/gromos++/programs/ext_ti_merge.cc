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
 * @file ext_ti_merge.cc
 * combine and reweight extended TI predictions
 */

/**
 * @page programs Program Documentation
 *
 * @anchor ext_ti_merge
 * @section ext_ti_merge combine and reweight extended TI predictions
 * @author @ref mp
 * @date 29. 4. 2017
 *
 * Program ext_ti_merge combines TI curves predicted by @ref ext_ti_ana from several simulations at different @f$\lambda@f$ points by a linear re-weighting scheme. @f$\frac{\partial H}{\partial \lambda} @f$ values predicted at a range of @f$\lambda_P@f$ between two simulated points @f$\lambda_{S1}@f$ and @f$\lambda_{S2}@f$ are reweighted in the following way:
 @f[ w_{S1} = \frac{\lambda_P-\lambda_{S2}}{\lambda_{S1}-\lambda_{S2}},  w_{S2} = \frac{\lambda_P-\lambda_{S1}}{\lambda_{S2}-\lambda_{S1}}@f]
 @f[\langle\frac{\partial H(r,t,\lambda_P)}{\partial \lambda} \rangle_{\lambda_P} =  w_{S1} \langle\frac{\partial H(r,t,\lambda_P)}{\partial \lambda} \rangle_{\lambda_{S1}}+ w_{S2} \langle\frac{\partial H(r,t,\lambda_P)}{\partial \lambda} \rangle_{\lambda_{S2}}@f]
 *
 * The integral (trapezoidal) of the final TI curve (and of its error values) is appended to the output file.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@files</td><td>&lt;data files (ext_ti_ana output)&gt; </td></tr>
 * <tr><td> [\@slam</td><td>&lt;lambdas of the simulation, optional if found in the header of the data files after '#SLAM '&gt;] </td></tr>
 * <tr><td> [\@lam_precision</td><td> &lt;lambda value precision in outfiles, default: 2&gt;] </td></tr>
 * <tr><td> [\@noerrors</td><td>&lt;do not read and use the error column which might be in the files&gt;] </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
  ext_ti_merge
    @files          predicted_l0.00.out predicted_l0.10.out predicte_l0.20.out
    @slam           0 0.1 0.2 
    @noerrors
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <algorithm> 

#include "../src/args/Arguments.h"
#include "../src/gromos/Exception.h"

using namespace args;
using namespace std;

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}
std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

// trim from start (in place)
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

double integrate(vector<double> x, vector<double> y);
double integrate_error(vector<double> x, vector<double> y);
vector<vector<double> > weight_function(ostream &os, vector<double> &slam, map<double, vector<double> > &vplam, map<double, vector<double> > &vdxdl, map<double, vector<double> > &vdxdl_err, bool no_errors);

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "files" << "slam" << "noerrors"  << "weightfunction" << "lam_precision";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@files            <data files (ext_ti_ana output)>\n";
  usage += "\t@noerrors         <do not read and use the error column which might be in the files>\n";
  usage += "\t[@slam            <lambdas of the simulation, optional if found in the header of the data files after '#SLAM '>]\n";
  usage += "\t[@lam_precision   <lambda value precision in outfiles, default: 2>]\n";
  //usage += "\t[@weightfunction      <only 'default' implemented so far>]\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    bool no_errors=false;
    if (args.count("noerrors") >= 0) no_errors=true;

    if (args.count("files") <= 1)
      throw gromos::Exception(argv[0], "At least two data files required!\n" +
            usage);
    
    string weightfunction_type = args.getValue<string>("weightfunction", false, "default");

    // lambdas given via flag?
    bool read_slam = false;
    vector<double> slam;
    if (args.count("slam") > 0) {
      if (args.count("slam") != args.count("files")) {
        throw gromos::Exception(argv[0],
              "@slam needs to provide one lambda value for each given file");
      } else {
        read_slam=true;
        for (Arguments::const_iterator iter = args.lower_bound("slam"),
            to = args.upper_bound("slam"); iter != to; ++iter) {
            double l = atof(iter->second.c_str());
            slam.push_back(l);
        }
      }
    }

    int lam_precision = args.getValue<int>("lam_precision",false, 2);

    map<double, vector<double> > vplam, vdxdl, vdxdl_err;
    
    // loop over the files and store data
    string line;
    stringstream linestream;    
    int fcnt=0;
    for (Arguments::const_iterator
      iter = args.lower_bound("files"),
            to = args.upper_bound("files") ;
            iter != to; ++iter, fcnt++) {

      ifstream file(iter->second.c_str());
      if (!file.good()) {
        throw gromos::Exception(argv[0], "Could not open file '" + iter->second + "'");
      }
      if (!file.is_open()) {
        throw gromos::Exception(argv[0], "Could not open file '" + iter->second + "'");
      }

      
      vector<double> plam, dxdl, dxdl_err;
      bool read_slam_file=false;
      do {
        getline(file, line, '\n');
        if (!file.eof()) {
          trim(line);
          linestream.clear();
          linestream.str(line);
          std::string match="#SLAM";
          double slami, plami, dxdli, dxdl_erri;
          if (!read_slam && !read_slam_file && line.find(match) != std::string::npos) {
            linestream >> match >> slami;
            slam.push_back(slami);
            read_slam_file=true;
          } else if (line[0] != '#') {
            linestream >> plami >> dxdli;
            
            // we couldn't read two doubles
            if (linestream.fail())
                throw gromos::Exception(argv[0], "can't read this file format, all lines not starting with a '#' should contain 2 or (with errors) 3 columns of floats");
            plam.push_back(plami);
            dxdl.push_back(dxdli);
                
            // is there also a value for the error?
            if (!linestream.eof() && !no_errors) {
              linestream >> dxdl_erri;
              dxdl_err.push_back(dxdl_erri);
            } else {
              no_errors=true;
            }
            if (linestream.fail()) { 
              no_errors=true;
            }
          }
        }
      } while (!file.eof());
        if (!read_slam_file && !read_slam)
            throw gromos::Exception(argv[0], "You have to provide the simulated lambda values either with @slam or in the file's header after '#SLAM '");

      vplam[slam[fcnt]]=plam;
      vdxdl[slam[fcnt]]=dxdl;
      vdxdl_err[slam[fcnt]]=dxdl_err;

      file.close();
    }

    // CALCULATIONS AND OUTPUT
    cout << "# Interpolating from data simulated at lambdas: ";
    for (vector<double>::const_iterator
      iter = slam.begin(),
            to = slam.end();
            iter != to; ++iter) {
      cout <<  *iter << " ";
    }
    cout << endl;
    
    std::vector<vector<double> > prediction;
    if (weightfunction_type=="default")
         prediction = weight_function(cout, slam, vplam, vdxdl, vdxdl_err, no_errors);
    else
         throw gromos::Exception(argv[0], "Unknown weightfunction "+weightfunction_type);
      
    
    double integral=integrate(prediction[0],prediction[1]);
    double interr=0;
    if (!no_errors) interr=integrate_error(prediction[0],prediction[2]);
    
    
    for (unsigned int i=0; i < prediction[0].size(); i++) {
      cout << setw(4)<< setprecision(lam_precision) << fixed<< prediction[0][i]<< " " << setprecision(6)<< fixed << prediction[1][i];
      if (!no_errors) cout << " " <<setprecision(6)<< fixed <<  prediction[2][i];
      cout << endl;
    }
    cout << "# integral: " << integral;
    if (!no_errors) cout << " +/- " << interr;
    cout << endl;



  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

double integrate(vector<double> x, vector<double> y) {
  double integral=0;
  
  if (x.size() != y.size()) throw gromos::Exception("integrate", "Vectors are not the same length!");
  
  for (unsigned int i=0; i<x.size()-1; i++) {
    integral+=0.5 * (y[i+1]+y[i])*(x[i+1]-x[i]);
  }
  return integral;
}

double integrate_error(vector<double> x, vector<double> y) {
  double integral=0;
  int N=x.size();
  if (N != y.size()) throw gromos::Exception("integrate", "Vectors are not the same length!");
  integral = 0.25 * (x[1]-x[0]) * (x[1]-x[0]) * y[0] * y[0];
  for (unsigned int i=1; i<N-2; i++) {
    integral+= y[i] * y[i] * (x[i+1]-x[i-1]) * (x[i+1]-x[i-1]) / 4;
  }
  integral += 0.25 * (x[N-1]-x[N-2]) * (x[N-1]-x[N-2]) * y[N-1] * y[N-1];
  return sqrt(integral);
}

vector<vector<double> > weight_function(ostream &os, vector<double> &slam, map<double, vector<double> > &vplam, map<double, vector<double> > &vdxdl, map<double, vector<double> > &vdxdl_err, bool no_errors) {
    
    sort(slam.begin(), slam.end());
    vector<vector<double> > pred_out(3);
    
    for (unsigned int i=0; i < slam.size()-1; i++) {
      double sl1=slam[i], sl2=slam[i+1];
      if (sl1==sl2)
        throw gromos::Exception("ext_ti_merge", "Two slam values are the same!\n");      
     
      // get the relevant predicted lambdas from i and i+1 in the range between i and i+1
      vector<vector<double> > plam_tmp(2), dxdl_tmp(2), dxdl_err_tmp(2);
      for (int s=0; s<=1; s++) {
        for (unsigned int j=0; j<vplam[slam[i+s]].size(); j++) {
          if (vplam[slam[i+s]][j]>=sl1 && vplam[slam[i+s]][j] <sl2) {
            plam_tmp[s].push_back(vplam[slam[i+s]][j]);
            dxdl_tmp[s].push_back(vdxdl[slam[i+s]][j]);
            if (!no_errors) dxdl_err_tmp[s].push_back(vdxdl_err[slam[i+s]][j]);
          }
        }        
      }
      if (plam_tmp[0] != plam_tmp[1]) {
        cerr << "ext_ti_merge: predicted lambdas between " << sl1 << " and " << sl2 << " are not the same, skipping.\n";
        continue;
      } else if (plam_tmp[0].size() == 0) {
        cerr << "ext_ti_merge: no predicted lambdas between " << sl1 << " and " << sl2 << ".\n";
        continue;
      } else {
        for (unsigned int j=0; j<plam_tmp[0].size(); j++) {
          double ws1 = (plam_tmp[0][j]-sl2)/(sl1-sl2);
          double ws2 = (plam_tmp[0][j]-sl1)/(sl2-sl1);
          double predict = dxdl_tmp[0][j]*ws1 + dxdl_tmp[1][j]*ws2;
          double error;
          if (!no_errors) {
            error = sqrt((dxdl_err_tmp[0][j]*dxdl_err_tmp[0][j]*ws1*ws1) + (dxdl_err_tmp[1][j]*dxdl_err_tmp[1][j]*ws2*ws2));
          }
          pred_out[0].push_back(plam_tmp[0][j]);
          pred_out[1].push_back(predict);
          if (!no_errors) pred_out[2].push_back(error);
        }
      }
    }
    
    // if the last simulated lambda is in the predicted ones (which it should be)
    // we still have to add it
    double last_slam=slam[slam.size()-1];
    std::vector<double>::iterator it=std::find(vplam[last_slam].begin(), vplam[last_slam].end(), last_slam);
    if (it != vplam[last_slam].end()) {
      int index=std::distance(vplam[last_slam].begin(), it);
      pred_out[0].push_back(vplam[last_slam][index]);
      pred_out[1].push_back(vdxdl[last_slam][index]);
      if (!no_errors) pred_out[2].push_back(vdxdl_err[last_slam][index]);
    }
    
    return pred_out;
}
