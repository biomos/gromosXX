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
 * @file bar.cc
 * Compute a free-energy difference using Bennett's acceptance ratio (BAR) method
 */

/**
 * @page programs Program Documentation
 *
 * @anchor bar
 * @section bar Compute a free-energy difference using Bennett's acceptance ratio (BAR) method
 * @author @ref co
 * @date 16.8.2017
 *
 * Program bar calculates free energy differences between (at least) two states using 
 * Bennett's Acceptance Ratio [Bennett, J. Comput. Phys. 22 (1976) 245 - 268]. The 
 * free energy between two states, i and j, is given by
 *
 * @f[ \Delta G(\lambda_i \rightarrow \lambda_j) 
 *     = k_B T \ln \frac{\langle f(E(\lambda_i) - E(\lambda_j) + C)\rangle_{\lambda_j}}
 *          {\langle f(E(\lambda_j)-E(\lambda_i) - C )\rangle_{\lambda_i}}  + C @f]
 *
 * where @f$ f(x) @f$ denotes the Fermi function
 *
 * @f[ f(x) = \frac{1}{1+exp(x/k_B T)} @f]
 * 
 * and
 *
 * @f[ C = k_B T ln \frac{N_j}{N_i} + \Delta G(\lambda_i \rightarrow \lambda_j). @f]
 *
 * These equations are solved self-consistently, using a numerically stable 
 * implementation adopted from [Shirts et. al. Phys. Rev. Lett. 91 (2003) 140601]. The
 * convergence criterion (relative change in free energy between iterations; 
 * \@convergence) and maximum number of iterations (\@maxiterations) can be specified.
 *
 * Time series of the energies (with lenght @f$ N_i @f$ and @f$ N_j @f$, respectively) 
 * are read in from one file per simulated state, which also contains the energies of 
 * the neighbouring states. These files may be generated using 
 * program @ref ext_ti_ana with option \@bar_data.
 *
 * Error estimates are standardly determined from the ensemble averages in the 
 * equations above. Optionally, a bootstrap error can be computed, where the 
 * calculation is repeated the indicated number of times (option \@bootstrap), with 
 * random samples of the original time series. The standard deviation of the bootstrap 
 * estimates is reported.
 * 
 * The program also computes the overlap integral from distributions of the energy
 * differences @f$ P_{i}(\Delta E) @f$ and @f$ P_{j}(\Delta E) @f$, using
 *
 * @f[ OI = 2 \sum_{\Delta E} \frac{P_{i}(\Delta E)P_{j}(\Delta E)}
 *          {P_{i}(\Delta E)+P_{j}(\Delta E)} @f] 
 *
 * with the sum running over all bins of the distribution. The distributions may be 
 * written out to separate files using option \@printdist.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@files</td><td>&lt;data files (@ref ext_ti_ana output)&gt; </td></tr>
 * <tr><td> \@temp</td><td>&lt;temperature&gt; </td></tr>
 * <tr><td> [\@maxiterations</td><td>&lt;maximum number of iterations (default: 500)&gt;] </td></tr>
 * <tr><td> [\@convergence</td><td>&lt;relative change in free energy for convergence (default: 1E-5)&gt;] </td></tr>
 * <tr><td> [\@printdist</td><td>&lt;write out distributions&gt;]</td></tr>
 * <tr><td> [\@bootstrap</td><td>&lt;<number of bootstrap estimates for error estimates (default: 0)&gt;] </td></tr>
 * </table>
 * 
 * Example:
 * @verbatim
  bar
    @files           bar_data_0.0.dat bar_data_0.20.dat bar_data_0.40.dat
    @temp            300
    @maxiterations   500
    @convergence     1E-11
    @printdist
    @bootstrap       100
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <algorithm> 

#include "../src/args/Arguments.h"
#include "../src/gmath/Stat.h"
#include "../src/gmath/Physics.h"
#include "../src/gmath/Distribution.h"
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

class barData
{
public:
  barData(){};
  
  double slam;
  std::vector<double> plam;
  std::vector<double> E_s;
  std::vector< std::vector<double> > E_p;

  barData(const barData &d){
    slam = d.slam;
    plam = d.plam;
    E_s = d.E_s;
    E_p = d.E_p;
  }
  
};

double do_bar(std::vector<double> &E_ii, std::vector<double> &E_ij,
	      std::vector<double> &E_jj, std::vector<double> &E_ji,
	      double & error,
	      int maxiter, double conv_eps);

double compute_overlap(std::vector<double> &E_ii, std::vector<double> &E_ij,
		       std::vector<double> &E_jj, std::vector<double> &E_ji,
		       int print);


  
int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "files" << "maxiterations" << "convergence" << "bootstrap" << "temp" << "printdist";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@files            <data files (ext_ti_ana bar_data output)>\n";
  usage += "\t@temp             <temperature>\n";
  usage += "\t[@maxiterations   <maximum number of iterations (default: 500)>]\n";
  usage += "\t[@convergence     <relative change in free energy for convergence (default: 1E-5)>]\n";
  usage += "\t[@printdist       <write out distributions>]\n";
  usage += "\t[@bootstrap       <number of bootstrap estimates for error estimates (default: 0)>]\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    int max_iter = args.getValue<int>("maxiterations", false, 500);
    double conv_eps = args.getValue<double>("convergence", false, 1E-5);
    double temp = args.getValue<double>("temp", true, 300);

    double kBT = gmath::physConst.get_boltzmann() * temp;
    bool printdist = false;
    if(args.count("printdist")>=0) printdist = true;
 
    int bootstrap = args.getValue<int>("bootstrap", false, 0);
    // initialize random number generator
    int ti=time(NULL);
	
    if (args.count("files") <= 1)
      throw gromos::Exception(argv[0], "At least two data files required!\n" +
            usage);

    
    std::vector<barData> bar_data;
    
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
      
      int nr_bar_lam = 0;
      double slam, plam, E_s, E_p;
      barData bd;

      do{ getline(file, line, '\n'); } while(line[0]=='#');
      trim(line);
      linestream.clear();
      linestream.str(line);
      linestream >> nr_bar_lam;
      if(linestream.fail())
	throw gromos::Exception(argv[0], "failed to read number of predicted lambda values");
      
      do{ getline(file, line, '\n'); } while(line[0]=='#');
      trim(line);
      linestream.clear();
      linestream.str(line);
      linestream >> bd.slam;
      for(int i=0; i< nr_bar_lam; i++) {
	linestream >> plam;
	bd.plam.push_back(plam);
      }
      if(linestream.fail())
	throw gromos::Exception(argv[0], "failed to read slam and or plam values");
      
      bd.E_p.resize(nr_bar_lam);
      do{ getline(file, line, '\n'); } while(line[0]=='#');
      while(!file.eof()){ 
	trim(line);
	linestream.clear();
	linestream.str(line);
	linestream >> E_s;
	//all energies are divided by kT
	bd.E_s.push_back(E_s / kBT);
	for(int i=0; i< nr_bar_lam; i++) {
	  linestream >> E_p;
	  bd.E_p[i].push_back(E_p / kBT);
	}
	if(linestream.fail() && !file.eof())
	  throw gromos::Exception(argv[0], "failed to read E_s and or E_p values");
	do{ getline(file, line, '\n'); } while(line[0]=='#');
      }
      for(int i=0; i< nr_bar_lam; i++){
	if(bd.E_s.size() != bd.E_p[i].size()){
	  throw gromos::Exception(argv[0], "stored different amount of E_s than E_p values");
	}
      }
      
      std::cerr << "read file " << iter->second << " with " << nr_bar_lam << " plam values and " << bd.E_s.size() << " data_points" << endl;
      
      bar_data.push_back(bd);
      
      file.close();
    }
    
    std::vector<double> lam_s;
    for(unsigned int i=0; i< bar_data.size(); i++){
      lam_s.push_back(bar_data[i].slam);
    }
    std::sort(lam_s.begin(), lam_s.end());

    //for(int i = 0; i < lam_s.size(); ++i)
    //  {
    //	std::cerr << "lam_s " << i << " : " << lam_s[i] << endl;
    // }
    
    double DG_tot = 0.0;
    double err_tot = 0.0;
    double boot_err_tot = 0.0;
    
    std::cout << "#" 
	      << setw(7) << "lam_i"
	      << setw(8) << "lam_j"
	      << setw(15) << "DG"
	      << setw(15) << "ee";
    if(bootstrap) std::cout << setw(15) << "bootstrap ee";
    std::cout << setw(15) << "overlap int."
	      << endl;

    int nr_intervals = bar_data.size()-1;
    for(int k=0; k< nr_intervals; k++){
      //std::cerr << "interval " << k << endl;
      
      //find matching indices to bar_data
      int i = -1;
      int j = -1;
      int i_j = -1;
      int j_i = -1;
      for(int l=0; l< bar_data.size(); l++){
	if(bar_data[l].slam == lam_s[k]) {
	  i = l;
	  for(int m =0; m < bar_data[l].plam.size(); m++){
	    if(bar_data[l].plam[m] == lam_s[k+1]) i_j = m;
	  }
	}
	if(bar_data[l].slam == lam_s[k+1]) {
	  j = l;
	  for(int m =0; m < bar_data[l].plam.size(); m++){
	    if(bar_data[l].plam[m] == lam_s[k]) j_i = m;
	  }
	}
      }
      if(i==-1 || j==-1)	
	throw gromos::Exception(argv[0], "could not find lambda_s value in data");
      if(i_j == -1){
	stringstream ss;
	ss << "Could not find data for plam = " << lam_s[k+1] << " in data file for slam = " << lam_s[k] << endl;
	throw gromos::Exception(argv[0], ss.str());
      }
      if(j_i == -1){
	stringstream ss;
	ss << "Could not find data for plam = " << lam_s[k] << " in data file for slam = " << lam_s[k+1] << endl;
	throw gromos::Exception(argv[0], ss.str());
      }
      //std::cerr << "lam_i " << lam_s[k] << " index " << i << endl;
      //std::cerr << "lam_j " << lam_s[k+1] << " index " << j << endl;
      //std::cerr << "lam_j in i " << lam_s[k+1] << " index " << i_j << endl;
      //std::cerr << "lam_i in j " << lam_s[k] << " index " << j_i << endl;

      double error=0.0;
      
      double DG = do_bar(bar_data[i].E_s, bar_data[i].E_p[i_j],
			 bar_data[j].E_s, bar_data[j].E_p[j_i],
			 error,
			 max_iter, conv_eps);

      //std::cerr << "result " << DG * kBT << " +/- " << error << std::endl;

      DG *= kBT;
      error *= kBT;
      
      DG_tot += DG;
      err_tot += error*error;

      double boot_std;
     
      if(bootstrap) {
	unsigned int seed = ti+k;
	
	int length_i = bar_data[i].E_s.size();
	int length_j = bar_data[j].E_s.size();

	gmath::Stat<double> final_boot;
	for(int ii=0;ii<bootstrap;ii++){
	  
	  std::vector<double> E_s_i, E_p_i, E_s_j, E_p_j;
	  
	  for(int jj=0;jj<length_i;jj++){
	    double r;
	    r=rand_r(&seed);
	    int index_i = int((r*length_i)/RAND_MAX);
	    	    
	    E_s_i.push_back(bar_data[i].E_s[index_i]);
	    E_p_i.push_back(bar_data[i].E_p[i_j][index_i]);
	  }
	  for(int jj=0; jj<length_j; jj++){
	    double r;
	    r=rand_r(&seed);
	    int index_j = int((r*length_j)/RAND_MAX);
	    
	    E_s_j.push_back(bar_data[j].E_s[index_j]);
	    E_p_j.push_back(bar_data[j].E_p[j_i][index_j]);
	  }
	  // don't compute the statistical error if bootstrapping
	  double err= -1.0;
	  double dg_boot = do_bar(E_s_i, E_p_i,
				  E_s_j, E_p_j,
				  err,
				  max_iter, conv_eps);
	  final_boot.addval(dg_boot*kBT);
	}
	boot_std = gmath::Stat<double>::covariance(final_boot, final_boot);
	boot_err_tot += boot_std;
      }
      
	
      // do we want to print the distributions?
      int print = -1;
      if(printdist) print = k;

      // calculate the overlap integral and print distributions if needed      
      double oi = compute_overlap(bar_data[i].E_s, bar_data[i].E_p[i_j],
				  bar_data[j].E_s, bar_data[j].E_p[j_i],
				  print);
      
      //std::cerr << "oi " << oi << endl;

      std::cout << setw(8) << fixed << setprecision(3) << lam_s[k]
		<< setw(8) << fixed << setprecision(3) << lam_s[k+1]
		<< setw(15) << setprecision(5) << DG
		<< " +/- "
		<< setw(10) << setprecision(5) << error;
      if(bootstrap) std::cout << setw(15) << setprecision(5) << sqrt(boot_std);
      std::cout << setw(15) << setprecision(5) << oi
		<< endl;
    }
    std::cout << endl
	      << "# total" 
	      << setw(24) << setprecision(5) << DG_tot
	      << " +/- "
	      << setw(10) << setprecision(5) << sqrt(err_tot);
    if(bootstrap) std::cout << setw(15) << setprecision(5) << sqrt(boot_err_tot);
    std::cout << endl;
    
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}



double do_bar(std::vector<double> &E_ii, std::vector<double> &E_ij,
	      std::vector<double> &E_jj, std::vector<double> &E_ji,
	      double &error,
	      int maxiter, double conv_eps){
  
  /*
    We will iterate to self consistency. The ensemble averages are computed in a 
    numerically stable way see [Shirts et. al. Phys. Rev. Lett. 91 (2003) 140601]
   */  
  int k=0;
  int N_i = E_ii.size();
  int N_j = E_jj.size();
  double M = log(double(N_i) / double(N_j));

  double reldiff = 1E15;
  double DG= 0;
  double DGprev= 0;
  
  // determine maximal forward energy difference and minimal backward energy difference
  // these won't change in the iterations so we don't have to determine them all the
  // time
  // For the forward guess, we will need the maximum of the forward energy difference
  // For the backward guess, we will need the maximum of the forward energy difference
  double max_ij = -1E15;
  double max_ji = -1E15;

  for(int ni = 0; ni < N_i; ni++)
    if((E_ij[ni] - E_ii[ni]) > max_ij) max_ij = E_ij[ni] - E_ii[ni];

  for(int nj = 0; nj < N_j; nj++)
    if((E_jj[nj] - E_ji[nj]) > max_ji) max_ji = E_jj[nj] - E_ji[nj];

  // calculate the exponentials in advance to keep the exp() out of the iteration and data loop
  std::vector<double> exp_ij;
  std::vector<double> exp_ji;
  for(int ni = 0; ni < N_i; ni++)
    exp_ij.push_back(exp(E_ij[ni] - E_ii[ni] - max_ij));
  for(int nj = 0; nj < N_j; nj++)
    exp_ji.push_back(exp(E_jj[nj] - E_ji[nj] - max_ji));

  // Go into the iterations
  while(k < maxiter && reldiff > conv_eps){
 
    /* 
       DG = ln<f(xi)>_i / ln<f(xj)>_j + C
       where C = DG - M
       and f(x) is the Fermi function f(x) = 1 / (1+exp(x))

       we need the following quantities. 
       - from simulation at i we need 

         ln<f(xi)>_i, with xi = E_ij - E_ii - DG + M

	 we will use a numerically stable routing to compute ln<exp(val)>, also used 
	 in the function lnexpave() of gmath::Stat. We write
	 ln<exp(log(f(xi))>, such that val = log(f(xi)). Define xi* as the 
         maximum of all xi

	 log(f(xi)) = log ( 1 / (1 + exp(xi)) )
                    = -log( 1 + exp(xi) )
                    = -log( exp(xi*) ( exp(-xi*) + exp(xi - xi*) ) )
                    = -xi* - log( exp(-xi*) + exp(xi - xi*) )

         filling in xi = E_ij - E_ii - DG + M and xi* = max_ij - DG + M = mx_ij
	 we get
	 
	 log(f(xi)) = -mx_ij - log(exp(-mx_ij) + exp(E_ij - E_ii - max_ij) )

	 or, using the previously computed exponentials, 

	 log(f(xi)) = -mx_ij -log(exp_mx_ij + exp_ij)

	 which we use to compute lnsum_ij, such that ln<f(xi)>_i = lnsum_ij - ln(N_i)
         the term ln(N_i) later cancels against the -M in C.

       - from simulation at j we need 

         ln<f(xj)>_j, with xj = E_ji - E_jj + DG - M 

	 we will use a numerically stable routing to compute ln<exp(val)>, also used 
	 in the function lnexpave() of gmath::Stat. We write
	 ln<exp(log(f(xj))>, such that val = log(f(xj))

	 We use f(xj) = f(-xj) exp(-xj) and define xj* as the maximum of all (-xj)
 
	 log(f(xj)) = log( f(-xj) ) - xj
	            = -log( 1 + exp(-xj) ) - xj
                    = -log( exp(xj*) ( exp(-xj*) + exp(-xj - xj*) ) ) - xj
                    = -xj* - log( exp(-xj*) + exp(-xj - xj*)) - xj

         filling in xj = E_ji - E_ii + DG - M and xj* = max_ji - DG + M = mx_ji
	 we get

	 log(f(xj)) = -mx_ji - log(exp(-mx_ji) + exp(E_jj - E_ji - max_ji) ) 
                                                          + E_jj - E_ji - DG + M

	 or, using the previously computed exponentials, and taking DG and M out of mx_ji

	 log(f(xj)) = -max_ji - log(exp_mx_ji + exp_ji) - E_jj - E_ji

	 which we use to compute lnsum_ji, such that ln<f(xj)>_j = lnsum_ji - ln(N_j)
         the term ln(N_j) later cancels against the -M in C.
    */

    DGprev = DG;

    double mx_ij = M + max_ij - DG;
    double exp_mx_ij = exp(-mx_ij);
    
    double mx_ji = M + max_ji - DG;
    double exp_mx_ji = exp(-mx_ji);

    double val=0;
    double lnsum_ij = -mx_ij - log(exp_mx_ij + exp_ij[0]);

    for(unsigned int ni = 1; ni < E_ii.size(); ++ni){
      val = -mx_ij - log(exp_mx_ij + exp_ij[ni]);
      lnsum_ij = std::max(lnsum_ij, val) + log(1.0 + exp(std::min(lnsum_ij, val) - std::max(lnsum_ij, val)));
    }

    double lnsum_ji = -max_ji - log(exp_mx_ji + exp_ji[0]) + E_jj[0] - E_ji[0];
    
    for(unsigned int nj = 1; nj < E_jj.size(); ++nj){
      val = -max_ji - log(exp_mx_ji + exp_ji[nj]) + E_jj[nj] - E_ji[nj];
      lnsum_ji = std::max(lnsum_ji, val) + log(1.0 + exp(std::min(lnsum_ji, val) - std::max(lnsum_ji, val)));
    }

    DG += lnsum_ji  - lnsum_ij;
    
    if(k>0) reldiff = abs((DG-DGprev)/DGprev);
    k++;
  }
  if (k==maxiter)
    std::cerr << "WARNING: Iteration finished because maximum number of iterations reached" << std::endl;

  
  // now for the error estimates, we do the same again
  // little trick: don't compute the error estimate if error is initialized negatively. 
  // Speeds up bootstrapping.
  if(!(error < 0)){

    // individual error estimates as in dfmult
    // in the iteration we computed the lnsum_ij ourselves to avoid copying the data back an forth.
    // here, we just store the appropriate 'val' in Delta_ij and Delta_ji, such that we can use the
    // features of the Stat class
    gmath::Stat<double> Delta_ij;
    gmath::Stat<double> Delta_ji;
  
    double mx_ij = M + max_ij - DG;
    double exp_mx_ij = exp(-mx_ij);
    
    double mx_ji = M + max_ji - DG;
    double exp_mx_ji = exp(-mx_ji);

    for(unsigned int ni = 0; ni < E_ii.size(); ++ni){
      Delta_ij.addval(-mx_ij - log(exp_mx_ij + exp_ij[ni]));
    }
    double d_ij = Delta_ij.lnexpave();
    
    // variance
    int sign_var_i = 1;
    double var_i = gmath::Stat<double>::lnexpcovariance(Delta_ij, Delta_ij, sign_var_i);
    if(sign_var_i < 0)
      throw gromos::Exception("bar", "Got negative variance!");
    // statistical inefficiency
    double si_i = gmath::Stat<double>::lnexp_stat_ineff(Delta_ij, Delta_ij);
    // error estimate contribution
    double d2i = var_i - log(double(N_i)) + si_i;
    double error_i2 = exp(d2i - 2*d_ij);
    
    for(unsigned int nj = 0; nj < E_jj.size(); ++nj){
      Delta_ji.addval(-max_ji - log(exp_mx_ji + exp_ji[nj]) + E_jj[nj] - E_ji[nj]);
    }
    double d_ji = Delta_ji.lnexpave();
    
    // variance
    int sign_var_j = 1;
    double var_j = gmath::Stat<double>::lnexpcovariance(Delta_ji, Delta_ji, sign_var_j);
    if(sign_var_j < 0)
      throw gromos::Exception("bar", "Got negative variance!");
    // statistical inefficiency
    double si_j = gmath::Stat<double>::lnexp_stat_ineff(Delta_ji, Delta_ji);
    // error estimate contribution
    double d2j = var_j - log(double(N_j)) + si_j;
    double error_j2 = exp(d2j - 2*d_ji);

    error = sqrt(error_i2 / (d_ij*d_ij) + error_j2 / (d_ji*d_ji));
  }
  
  return DG;
   
}

double compute_overlap(std::vector<double> &E_ii, std::vector<double> &E_ij,
		       std::vector<double> &E_jj, std::vector<double> &E_ji,
		       int print){

  gmath::Stat<double> Delta_ij;
  gmath::Stat<double> Delta_ji;
  
  for(unsigned int ni = 0; ni < E_ii.size(); ++ni){
    Delta_ij.addval(E_ij[ni] - E_ii[ni]);
  }
  for(unsigned int nj = 0; nj < E_jj.size(); ++nj){
    Delta_ji.addval(E_jj[nj] - E_ji[nj]);
  }
  
  double min, max;
  
  if(Delta_ij.min() < Delta_ji.min()) min = Delta_ij.min();
  else min = Delta_ji.min();
  if(Delta_ij.max() > Delta_ji.max()) max = Delta_ij.max();
  else max = Delta_ji.max();
  int nbins = 100;
  
  Delta_ij.dist_init(min, max, nbins, false);
  Delta_ji.dist_init(min, max, nbins, false);
  
  double oi=0;
  
  for(int n=0; n< nbins; n++){
    double Pi = Delta_ij.distribution()[n]/double(Delta_ij.distribution().nVal());
    double Pj = Delta_ji.distribution()[n]/double(Delta_ji.distribution().nVal());
    if(Pi+Pj) oi += (Pi*Pj)/(Pi+Pj);
  }

  if(print >=0){
    int k=print;
    stringstream ss;
    ss << "dist_" << k << "_" << k+1 << ".dat";
    ofstream fout(ss.str().c_str());
    fout << "# Distribution from state " << k << " to state " << k+1 << endl;
    fout << "#" << endl;
    fout << "# DeltaE    P(DeltaE)" << endl;
    Delta_ij.distribution().write_normalized(fout);
    fout << "# Distribution from state " << k+1 << " to state " << k << endl;
    fout << "#" << endl;
    fout << "# DeltaE    P(-DeltaE)" << endl;
    Delta_ji.distribution().write_normalized(fout);
    fout.close();
  }
  
  return 2*oi;
}

