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
 * @file contactnum.cc
 * calculates atom-positional root-mean-square deviations
 */
/**
 * @page programs Program Documentation
 *
 * @anchor contactnum
 * @section contactnum calculates the number of contacts (distances smaller than a given cutoff) between two groups of atoms
 * @author @ref mp
 * @date 14.1.13
 *
 * The program calculates a descriptor of the the number of contacts between two 
 * sets of atoms according to the following equation:
 * ..
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gathermethod&gt;] </td></tr>
 * <tr><td> \@time</td><td>&lt;@ref utils::Time "time and dt"&gt; </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier two atomspecifiers separated by whitespace&gt; </td></tr>
 * <tr><td> [\@cutoffs</td><td>&lt;cutoff distances&gt;] </td></tr>
 * <tr><td> [\@exp</td><td>&lt;n and m exponents&gt;] </td></tr>
 * <tr><td> [\@printpairs</td><td>&lt;print pairs to file&gt;] </td></tr>
 * <tr><td> [\@excludemass</td><td>&lt;exclude atoms with this mass&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  ditrans
    @topo       ex.top
    @pbc        r
    @time       0 0.1
    @atoms  1:CA 2:CA
    @cutoffs  0.4 1.0 1.4
    @exp        2 24
    @excludemass 1.008
    @printpairs pairs.txt
    @traj       ex.trc
    @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/args/ReferenceParser.h"
#include "../src/fit/Reference.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
#include "../src/gmath/Vec.h"
#include "../src/gromos/Exception.h"


using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;
using namespace std;
using namespace gmath;

double switchingfunction(double rdist,int nn,int mm);
double fastpow(double base, int exp);
  
int main(int argc, char **argv){
  Argument_List knowns; 
  knowns << "topo" << "traj" << "atoms" << "cutoffs" << "pbc" << "exp"
         << "time" << "printpairs" << "excludemass";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo         <molecular topology file>\n";
  usage += "\t@pbc          <boundary type> [<gathermethod>]\n";
  usage += "\t@time         <time and dt>\n";
  usage += "\t@atoms        <two groups of atoms>\n";
  usage += "\t@cutoffs      <cutoff distances >\n";
  usage += "\t[@exp         <n and m exponents, default: 2 24>]\n";
  usage += "\t[@printpairs  <print pairs to file>]\n";
  usage += "\t[@excludemass <exclude atoms with this mass/these masses>]\n";
  usage += "\t@traj         <trajectory files>\n";


  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);
  

  try{
    Arguments args(argc, argv, knowns, usage);

    // get simulation time
    Time time(args);
    
    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    
    if (args.count("pbc") <=0) {
      cout << "# WARNING: no @pbc argument found - no gathering will be done!";
    }
    
    bool writepairs=false;
    std::ofstream pairfile;
    if(args.count("printpairs")>=0) { 
      writepairs=true;
      std::string printpairs=args.getValue<string>("printpairs",false,"pairs.txt");
      pairfile.open(printpairs.c_str());
      pairfile << "# atom1 atom2 dist value\n";
    }
     
    vector<double> excludemass;
    if (args.count("excludemass") > 0) {
       Arguments::const_iterator iter = args.lower_bound("excludemass");
       Arguments::const_iterator to = args.upper_bound("excludemass");
       for (;iter!=to;iter++) {
        excludemass.push_back(atof(iter->second.c_str())); 
       }
    }
    else {
      excludemass.push_back(-1);
    }
    
     
    std::vector<int> exp = args.getValues<int>("exp", 2, false,
            Arguments::Default<int>() << 2 << 24);
    int nn= exp[0];
    int mm= exp[1];

    std::vector<double> cutoffs;
    if (args.count("cutoffs")>0) {
       Arguments::const_iterator iter = args.lower_bound("cutoffs");
       Arguments::const_iterator to = args.upper_bound("cutoffs");
       for(;iter!=to;iter++){
         cutoffs.push_back(atof(iter->second.c_str()));
       }
    }
    else {
      throw gromos::Exception("contactnum", "You need to give at least one cutoff!");
    }
    

    // read 2 atom groups
    AtomSpecifier atomgroup1(sys), atomgroup2(sys);
    if (args.count("atoms")>1) {
      {
       Arguments::const_iterator iter = args.lower_bound("atoms");
       Arguments::const_iterator to = args.upper_bound("atoms");
       
       atomgroup1.addSpecifier(iter->second.c_str());
       iter++;
       atomgroup2.addSpecifier(iter->second.c_str());
       iter++;
       if (iter!=to)
         throw gromos::Exception("contactnum", "Too many atom groups! Give two atom selections separated by whitespace!");
      }  
    }
    else {
      throw gromos::Exception("contactnum", "Too few atom groups! Give two atom selections separated by whitespace!");
    }
    
    // exclude atoms with given mass
    for (int i=atomgroup1.size()-1; i>=0 ; i--) {
      for(unsigned int j=0; j<excludemass.size(); j++) {
        if (atomgroup1.mass(i)== excludemass[j]) {
          //cerr << "atom " << atomgroup1.atom(i) << " mass " << atomgroup1.mass(i)<<endl;
          atomgroup1.removeAtom(i);
        }
      }
    }
    for (int i=atomgroup2.size()-1; i>=0 ; i--) {
      for(unsigned int j=0; j<excludemass.size(); j++) {
        if (atomgroup2.mass(i)== excludemass[j]){
          //cerr << "atom " << atomgroup2.atom(i) << " mass " << atomgroup2.mass(i) << endl;
          atomgroup2.removeAtom(i);
        }
      }
    }  
    if (atomgroup1.size()==0 || atomgroup2.size()==0)
        throw gromos::Exception("contactnum", "No atoms in one or both of the specified atom groups!");    
    
    cerr << "# " << atomgroup1.size() << " atoms in group 1" << endl;
    cerr << "# " << atomgroup2.size() << " atoms in group 2" << endl;
    cout  <<  "# time         ";
    // Parse boundary conditions for sys
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    
	for (unsigned int cut=0; cut<cutoffs.size(); cut++) {
	   cout <<setw(15) << cutoffs[cut];
	}
	cout << endl;
    
    // loop over all trajectories
    int numFrames = 0;
    // average contact number
    vector<double> ct_ave(cutoffs.size(),0.0);
    InG96 ic;
    for (Arguments::const_iterator iter=args.lower_bound("traj");
	                      iter!=args.upper_bound("traj"); ++iter){
      ic.open(iter->second);
      
      // loop over all frames
      while (!ic.eof()) {
          
        numFrames++;
        ic.select("ALL");
	    ic >> sys >> time;
        if(!sys.hasPos)
          throw gromos::Exception("contactnum",
                             "Unable to read POSITION(RED) block from "
                              "trajectory file.");
	
        cout.precision(2);
        cout << time;
	    for (unsigned int cut=0; cut<cutoffs.size(); cut++) {
	      double ct=0;
          for (unsigned int i=0; i< atomgroup1.size(); i++) {
            for (unsigned int j=0; j < atomgroup2.size(); j++) {
              double func;
              gmath::Vec v = atomgroup1.pos(i) - pbc->nearestImage(atomgroup1.pos(i), atomgroup2.pos(j), sys.box());
              double rdist=v.abs();
              func=switchingfunction(rdist/cutoffs[cut], nn, mm);
              ct+=func;
              if (writepairs) {
                pairfile << atomgroup1.mol(i)+1 <<":"<<atomgroup1.atom(i)+1<< " "<< atomgroup2.mol(j)+1 <<":"<<atomgroup2.atom(j)+1 << " " << rdist << " " << func <<" " <<ct <<endl;
              }
            }
          }  
          cout.precision(9);
          cout << setw(15) << ct ;
          ct_ave[cut] += ct;
        }
        cout << endl;
      }
      ic.close();
    }  
    cout <<"# averages     ";
    for (unsigned int i=0; i<cutoffs.size(); i++) {
      cout << setw(15) << ct_ave[i]/numFrames; 
    }
    cout << endl;
    if (writepairs) {  
      pairfile.close();
    }
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


inline
double fastpow(double base, int exp)
{
    if(exp<0){
      exp=-exp;
      base=1.0/base;
    }
    double result = 1.0;
    while (exp)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }

    return result;
}

// copied from Plumed-2.2.0 switching function "rational"
double switchingfunction(double rdist,int nn,int mm) {
      // Very small non-zero number
      const double epsilon(std::numeric_limits<double>::epsilon());
      
      double result;
      if(2*nn==mm){
// if 2*N==M, then (1.0-rdist^N)/(1.0-rdist^M) = 1.0/(1.0+rdist^N)
        double rNdist=fastpow(rdist,nn-1);
        double iden=1.0/(1+rNdist*rdist);
        result = iden;
      } else {
        if(rdist>(1.-100.0*epsilon) && rdist<(1+100.0*epsilon)){
           result=nn/mm;
        }else{
           double rNdist=fastpow(rdist,nn-1);
           double rMdist=fastpow(rdist,mm-1);
           double num = 1.-rNdist*rdist;
           double iden = 1./(1.-rMdist*rdist);
           double func = num*iden;
           result = func;
        }
      }
    return result;
}

