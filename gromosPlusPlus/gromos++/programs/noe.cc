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
 * @file noe.cc
 * Analysis of NOE distances over a trajectory
 */

/**
 * @page programs Program Documentation
 *
 * @anchor noe
 * @section noe Analysis of NOE distances over a trajectory
 * @author @ref mk
 * @date 15.8.2006
 *
 * Program noe calculates and averages atom-atom restraint distances for 
 * specified NOE distances over a molecular trajectory. The NOE distances are 
 * to be specified in a NOE specification file, that can be prepared with e.g. 
 * program @ref prep_noe.
 *
 * Program NOE will calculate the average distance according to 
 * @f$<r^{-p}>^{-1/p}@f$ for values of p=1, 3, 6. It will also calculate the 
 * deviations of these distances from the specified reference distances, 
 * @f$r_0@f$. This violations can be written to a time series file. The average
 * violation is calculated as the sum of positive  violations
 * (i.e. if @f$(<r^{-p}>^{-1/p} - r_0) > 0@f$) divided by the total number of
 * NOE distances considered in the analysis. The output of the program can be
 * further analysed using program @ref post_noe "post_noe".
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> \@noe</td><td>&lt;NOE specification file&gt; </td></tr>
 * <tr><td> \@time</td><td>&lt;@ref utils::Time "time and dt"&gt; (if absent, no time series is written)</td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  noe
    @topo   ex.top
    @pbc    r
    @noe    noe.spec
    @traj   ex.tr
 @endverbatim
 *
 * <hr>
 */
#include <cassert>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>

#include "../src/args/Arguments.h"
#include "../src/gio/Ginstream.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/StringTokenizer.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/utils/Noe.h"
#include "../src/gmath/Stat.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/groTime.h"
#include "../src/utils/Value.h"
#include "../src/gcore/System.h"
#include "../src/bound/Boundary.h"
#include "../src/gromos/Exception.h"
#include "../src/utils/groTime.h"

using namespace bound;
using namespace gcore;
using namespace args;
using namespace gio;
using namespace utils;
using namespace std;


int main(int argc,char *argv[]){


  // Usage string
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t@pbc    <boundary type> [<gather method>]\n";
  usage += "\t@noe    <NOE specification file>\n"; 
  usage += "\t[@time   <time and dt>]\n";
  usage += "\t@traj   <trajectory files>\n";

  // known arguments...
  Argument_List knowns;
  knowns << "topo" << "noe" << "pbc" << "time" << "traj";
    
  // prepare cout for formatted output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(3);
    
  try{

    // Getting arguments and checking if everything is known.
    Arguments args(argc,argv,knowns,usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    System refSys(it.system());

    // get simulation time
    Time time(args);

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    //parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);
    // Read in and create the NOE list
    Ginstream nf(args["noe"]);
    vector<string> buffer;
    nf.getblock(buffer);
    
    if(buffer[0]!="NOECALCSPEC")
      throw gromos::Exception("main","NOE file does not contain an NOECALCSPEC block!");
    if(buffer[buffer.size()-1].find("END")!=0)
      throw gromos::Exception("noe", "NOE file " + nf.name() +
				" is corrupted. No END in "+buffer[0]+
				" block. Got\n"
				+ buffer[buffer.size()-1]);

    // in noe all noes will be stored.
    vector<Noe*> noe;

    string line;
    StringTokenizer tok(buffer[1]);
    vector<string> dishdisc = tok.tokenize();
    double dish = atof(dishdisc[0].c_str());
    double disc = atof(dishdisc[1].c_str());
    for(unsigned int j=2; j< buffer.size()-1; j++){      
      noe.push_back(new Noe(sys, buffer[j], dish, disc));
    }
    
    // vectors to contain the r^-3 and r^-6 averages, rmsds and errors
    vector<vector<gmath::Stat<double> > > s, s3, s6;
    vector<vector<double> > av, av3, av6;
    vector<vector<double> > ee, ee3, ee6;
    vector<vector<double> > rmsd, rmsd3, rmsd6;
    
    // initialisation of storage space
    unsigned int num_noes = noe.size();
    s.resize(num_noes);
    s3.resize(num_noes);
    s6.resize(num_noes);   
    av.resize(num_noes);
    av3.resize(num_noes);
    av6.resize(num_noes);
    ee.resize(num_noes);
    ee3.resize(num_noes);
    ee6.resize(num_noes);
    rmsd.resize(num_noes);
    rmsd3.resize(num_noes);
    rmsd6.resize(num_noes);
    
    for(int i=0; i<int(noe.size()); ++i){
      int n = noe[i]->numDistances();
      s[i].resize(n);
      s3[i].resize(n);
      s6[i].resize(n);
      av[i].resize(n);
      av3[i].resize(n);
      av6[i].resize(n);
      ee[i].resize(n);
      ee3[i].resize(n);
      ee6[i].resize(n);
      rmsd[i].resize(n);
      rmsd3[i].resize(n);
      rmsd6[i].resize(n);
    }

    //spit out title block
    cout << "TITLE" << endl;
    cout << "NOE analysis according to: " << args["noe"] << endl;
    cout << nf.title();
    cout << "END" << endl;
 
    nf.close();
    
    //  open timeseries files if requested        
    ofstream timeseries;
    if (time.doSeries()) {
      timeseries.open("noets.out");
      timeseries.setf(ios::right, ios::adjustfield);
      timeseries.setf(ios::fixed, ios::floatfield);
    }
    
    // define input coordinate
    InG96 ic;
    
    int numFrames=0;
    
// loop over all trajectories
    for(Arguments::const_iterator 
	  iter=args.lower_bound("traj"),
	  to=args.upper_bound("traj");
	iter!=to; ++iter){
      
      // open file
      ic.open((iter->second).c_str());
      
      // loop over single trajectory
      while(!ic.eof()){
	numFrames++;
	ic >> sys;
        if (time.doSeries())
          ic >> time;
	(*pbc.*gathmethod)();
	
	// loop over noes
	for(int nr = 1, i=0; i < int(num_noes); ++i) {
          for (int ii = 0; ii < noe[i]->numDistances(); ++ii, ++nr) {
            // calculate distance and averages...
	    double distance=noe[i]->distance(ii);
            double idist3 = 1 / (distance*distance*distance);
            double idist6 = idist3 * idist3;
            s[i][ii].addval(distance);
            s3[i][ii].addval(idist3);
            s6[i][ii].addval(idist6);
            
            if (time.doSeries()) {
              // only do this when time series is asked
              // calculation of errors at every step is too slow.
              av[i][ii] += distance;
              av3[i][ii] += idist3;
              av6[i][ii] += idist6;
              double ave = av[i][ii]/numFrames;
              double ave3 = pow(av3[i][ii]/numFrames, -1.0/3.0);
              double ave6 = pow(av6[i][ii]/numFrames, -1.0/6.0);

              timeseries.precision(2);
              timeseries << setw(10) << time << setw(8) << nr;
              timeseries.precision(3);
              timeseries << setw(8) << ave - noe[i]->reference(ii)
                         << setw(8) << ave3 - noe[i]->reference(ii)
                         << setw(8) << ave6 - noe[i]->reference(ii)
                         << std::endl;
            }
	  }
        }
      }
      ic.close();
    }
    
    if (time.doSeries())
      timeseries.close();

    // calculate the averages
    for(int i=0; i < int(noe.size());  ++i)
      for(int ii=0; ii < noe[i]->numDistances(); ++ii){
        av[i][ii] = s[i][ii].ave();
        ee[i][ii] = s[i][ii].ee();
        rmsd[i][ii] = s[i][ii].rmsd();
        double ave3 = s3[i][ii].ave();
        av3[i][ii] = pow(ave3,-1.0/3.0);
        ee3[i][ii] = std::abs(pow(ave3, -4.0/3.0) / 3.0) * s3[i][ii].ee();
        rmsd3[i][ii] = std::abs(pow(ave3, -4.0/3.0) / 3.0) * s3[i][ii].rmsd();
        double ave6 = s6[i][ii].ave();
        av6[i][ii] = pow(ave6,-1.0/6.0);
        ee6[i][ii] = std::abs(pow(ave6, -7.0/6.0) / 6.0) * s6[i][ii].ee();
        rmsd6[i][ii] = std::abs(pow(ave6, -7.0/6.0) / 6.0) * s6[i][ii].rmsd();
      }
    

    // output the averages
    cout << "AVERAGE NOE\n";
    for(int i=0, nr=1; i < int(num_noes); ++i)
      for(int j=0; j < noe[i]->numDistances(); ++j, ++nr)
	cout << "# " << setw(4) << nr << " " << noe[i]->info(j) << endl;
    
    cout <<'#' << ' ' 
	 << setw(4) << "Nr."
	 << setw(10) << "av"
         << setw(10) << "av3"
         << setw(10) << "av6"
         << setw(10) << "rmsd"
         << setw(10) << "rmsd3"
         << setw(10) << "rmsd6" 
         << setw(10) << "ee"
         << setw(10) << "ee3"
         << setw(10) << "ee6" << endl;

    for(unsigned int i=0, nr=1; i<num_noes; ++i)
      for(unsigned int ii=0, nnum=noe[i]->numDistances(); ii<nnum;++ii, ++nr)
	cout << setw(6) << nr
	     << setw(10) << av[i][ii]
	     << setw(10) << av3[i][ii]
	     << setw(10) << av6[i][ii]
 	     << setw(10) << rmsd[i][ii]
	     << setw(10) << rmsd3[i][ii]
	     << setw(10) << rmsd6[i][ii]
	     << setw(10) << ee[i][ii]
	     << setw(10) << ee3[i][ii]
	     << setw(10) << ee6[i][ii]
	     << endl;

    
    // Now treat the violations
    
    // order the distances first
    vector<vector<int> > order(num_noes);
    for(unsigned int i=0; i<num_noes; ++i){
      order[i].push_back(0);
      for(unsigned int ii=1, nnum=noe[i]->numDistances(); ii<nnum;++ii)
	for(vector<int>::iterator
	      it=order[i].begin(),
	      to=order[i].end()+1;
	    it!=to;++it)
	  if(it==to-1)order[i].insert(it,ii);
	  else if(av3[i][ii]<av3[i][*it]){
	    order[i].insert(it, ii);
	    break;
	  }
    }

    double cdaver=0.0, daver=0.0, avresvio=0.0, avresvio3=0.0, avresvio6=0.0; int c=0; 
    // now output the NOE violations
    cout << "END\nNOE VIOLATIONS\n";
    for(unsigned int i=0, nr=1; i<num_noes; ++i)
      for(unsigned int ii=0, nnum=noe[i]->numReferences(); ii<nnum;++ii, ++nr)
	cout << "# " << setw(4) << nr << " " << noe[i]->info(ii) << endl;
    cout << "#\n";
    cout << "# r0 = reference length according to specification file\n";
    cout << "# " 
	 << setw(4) << "Nr."
	 << setw(10) << "r0"
	 << setw(10) << "<r> - r0"
         << setw(20) << "<r^-3>^-1/3 - r0"
         << setw(20) << "<r^-6>^-1/6 - r0" << endl;
    
    for(unsigned int i=0, nr=1; i<num_noes; ++i){
      for(unsigned int ii=0, nnum=noe[i]->numReferences(); ii<nnum;++ii, ++nr){
	double cd;
	
	cd=noe[i]->reference(ii);
	
	//averaging stuff
        daver+=	noe[i]->reference(ii);
        cdaver+=cd;
        c++;
        if (av[i][order[i][ii]] - cd > 0.0){
          avresvio+=av[i][order[i][ii]] - cd;
	}
        if (av3[i][order[i][ii]] - cd > 0.0){
          avresvio3+=av3[i][order[i][ii]] - cd;
	}
        if (av6[i][order[i][ii]] - cd > 0.0){
          avresvio6+=av6[i][order[i][ii]] - cd;
	}

	//the real printout
	cout <<	setw(6) << nr 
	     << setw(10) << cd 
	     << setw(10) << av[i][order[i][ii]] - cd
	     << setw(20) << av3[i][order[i][ii]] - cd 
	     << setw(20) << av6[i][order[i][ii]] - cd << endl;
      }
    }

    cout << "END\n";
    cout << "# AVERAGE r0: " << cdaver/c << endl;
    cout << "# AVERAGE RESTRAINT VIOLATION (<av-r0>): " << avresvio/c << endl;
    cout << "# AVERAGE RESTRAINT VIOLATION (<av3-r0>): " << avresvio3/c << endl;
    cout << "# AVERAGE RESTRAINT VIOLATION (<av6-r0>): " << avresvio6/c << endl;
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
