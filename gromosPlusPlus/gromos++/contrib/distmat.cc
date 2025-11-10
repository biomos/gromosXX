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
 * @file distmat.cc
 * Calculate a distance matrix between structures 
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor distmat
 * @section distmat Calculate a distance matrix between structures
 * @author @ref co
 * @date 12-07-07
 *
 * Program distmat calculates a similarity matrix for all pairs of structures
 * in a given trajectory file, very much like the program rmsdmat. Rather than
 * the rmsd, this program uses a quantity that only depends on intramolecular
 * distances. It is therefore independent of the rotational fit. The distance
 * between two structures, 1 and 2 is calculated as 
 * @f[ D = \sqrt{\frac{1}{N_{pairs}}\sum_{i,j} [r_{ij}(1) - r_{ij}(2)]^2} @f]
 * where the sum goes over all possible pairs within the selected set of atoms.
 * This matrix of can subsequently be used by program @ref cluster to perform a
 * conformational clustering. The matrix can be written out in human readable
 * form, or -to save disk space- in binary format. For efficiency reasons, the
 * D values are written in an integer format. The user can specify the
 * required precision of the RMSD values that are stored, if the precision is
 * less or equal to 4, the values are stored as unsigned short int, otherwise
 * as unsigned int.
 *
 * A selection of structures in the trajectory file to consider can be made 
 * using the options skip and stride. Specifying a reference structure allows 
 * the program @ref cluster to perform a forced clustering as well, requiring 
 * that the first cluster contains the reference structure, regardless of the 
 * cluster size.
 *
 * <B>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary conditions&gt; &lt;gather type&gt; </td></tr>
 * <tr><td> \@atomsdist</td><td>&lt;@ref AtomSpecifier "atoms" to consider&gt; </td></tr>
 * <tr><td> [\@skip</td><td>&lt;skip frames at beginning&gt;] </td></tr>
 * <tr><td> [\@stride</td><td>&lt;use only every step frame&gt;] </td></tr>
 * <tr><td> [\@human</td><td>(write the matrix in human readable form)] </td></tr>
 * <tr><td> [\@precision</td><td>&lt;number of digits in the matrix (default 4)&gt;] </td></tr>
 * <tr><td> [\@big</td><td>(when clustering more than ~50'000 structures)] </td></tr>
 * <tr><td> [\@ref</td><td>&lt;reference coordinates&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  rmsdmat
    @topo         ex.top
    @pbc          r
    @atomsdist    1:CA
    @skip         5
    @stride       10
    @human
    @precision    4
    @big
    @ref          exref.coo
    @traj         ex.tr
 @endverbatim
 *
 * <hr>
 */


#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/Boundary.h"
#include "../src/gcore/System.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Matrix.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/fit/AtomDistances.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gmath;
using namespace gcore;
using namespace bound;
using namespace gio;
using namespace utils;
using namespace args;
using namespace fit;

int main(int argc, char **argv){

  Argument_List knowns; 
  knowns << "topo" << "traj" << "pbc" << "ref" << "atomsdist"
         << "skip" << "stride" << "human" << "precision" << "big"; 

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo         <molecular topology file>\n";
  usage += "\t@pbc          <boundary conditions> <gather type>\n";
  usage += "\t@atomsdist    <atoms to consider for fit>\n";
  usage += "\t[@skip        <skip frames at beginning>]\n";
  usage += "\t[@stride      <use only every step frame>]\n";
  usage += "\t[@human       (write the matrix in human readable form)]\n";
  usage += "\t[@precision   <number of digits in the matrix (default 4)>]\n";
  usage += "\t[@big         (when clustering more than ~50'000 structures)]\n";
  usage += "\t[@ref         <reference coordinates>]\n";
  usage += "\t@traj         <trajectory files>\n";

  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);

  try{
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    

    // read the dist atoms
    AtomSpecifier distatoms(sys);
    {
      Arguments::const_iterator iter=args.lower_bound("atomsdist"),
	to=args.upper_bound("atomsdist");
      for( ; iter!=to; ++iter){
	distatoms.addSpecifier(iter->second);
      }
      if(distatoms.size()==0)
	throw gromos::Exception("rmsdmat",
				"No fit atoms given\n");
    }
    
    // for which atoms do we want to keep the coordinates
    AtomSpecifier atoms(distatoms);
    
    // if fitatoms != rmsdatoms keep lists of what to do
    vector<bool> dist_spec;
    // this is empty, because we don't distinguish between fit atoms and 
    // rmsd atoms. That means that we could have kicked it from the class?
    AtomDistances ad(dist_spec);
    
    int skip = args.getValue<int>("skip", false, 0);
    int stride= args.getValue<int>("stride", false, 1);

    // read the precision
    int precision=int(pow(10.0, args.getValue<double>("precision", false, 4.0)));
    
    System refSys(it.system());

    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // GatherParser
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

    // create the vector to store the trajectory
    vector< vector < Vec > > traj;
    vector< Vec > frame(atoms.size());
    
    // read reference coordinates...
    InG96 ic;
    if(args.count("ref") > 0){
      ic.open(args["ref"]);
    }
    else{
      ic.open(args.lower_bound("traj")->second);
    }
    ic >> sys;
    ic.close();

    (*pbc.*gathmethod)();

    // put it in the trajectory
    for(int i=0; i < atoms.size(); ++i){
      frame[i] = atoms.pos(i) ;
    }
    traj.push_back(frame);
    
    int framenum=0;
    
    // loop over all trajectories
    for(Arguments::const_iterator iter=args.lower_bound("traj");
      iter!=args.upper_bound("traj"); ++iter){

      ic.open(iter->second);

      // loop over all frames
      while(!ic.eof()){
	ic >> sys;
	
	if (!((framenum - skip) % stride)){
	  
	  //pbc call
	  (*pbc.*gathmethod)();

	  for(int i=0; i < atoms.size(); ++i){
	    frame[i] = atoms.pos(i);
	  }
	  
          // store coordinates from sys in traj
	  traj.push_back(frame);
        }
	framenum++;
      }
      ic.close();
    }

    // everything is in the thing; create the thingy
    // open a file
    ofstream fout;
    bool human=false;
    if(args.count("human") >= 0){
      fout.open("DISTMAT.dat");
      human=true;
    }
    else{
      fout.open("DISTMAT.bin", ios::out | ios::binary);
    }
    
    // make a double loop
    int num = traj.size();
    Matrix rot(3,3,0);
    double dist = 0;
    Matrix unit(3,3,0);
    for(size_t i=0; i<3; i++) unit(i,i)=1;

    cout << "Read " << num << " out of " << framenum 
	 << " structures from trajectory" << endl;
    
    if(human){
      fout << "TITLE\n"
	   << "\trmsd-matrix for " << num-1 << " + 1 (ref) = " 
	   << num << " structures\n"
	   << "END\n"
	   << "RMSDMAT\n"
	   << "# number of frames   skip   stride\n"
	   << num << "\t" << skip << "\t" << stride << "\n"
	   << "# precision\n"
	   << precision << "\n";
      

      for(int i=0; i< num; ++i){
	for(int j=i+1; j < num; ++j){

	  dist = ad.dist(traj[i], traj[j]);
	  fout << setw(8) << i 
	       << setw(8) << j 
	       << setw(8) << unsigned(dist * precision)
	       << endl;
	}
      }
      fout << "END\n";
    }
    else{
      fout.write((char*)&num, sizeof(int));
      fout.write((char*)&skip, sizeof(int));
      fout.write((char*)&stride, sizeof(int));
      fout.write((char*)&precision, sizeof(int));
      
      if(precision < 1e5 && args.count("big") == -1){
	std::cout << "using 'unsigned short' as format" << std::endl;

	typedef unsigned short ushort;
	
	ushort irmsd;
	for(int i=0; i< num; ++i){
	  for(int j=i+1; j < num; ++j){
	    dist = ad.dist(traj[i], traj[j]);
	    irmsd=ushort(dist*precision);
	    
	    fout.write((char*)&irmsd, sizeof(ushort));
	  }
	}
      }
      else{
	std::cout << "using 'unsigned int' as format" << std::endl;
	unsigned irmsd;
	for(int i=0; i< num; ++i){
	  for(int j=i+1; j < num; ++j){
	    
	    dist = ad.dist(traj[i], traj[j]);
	    irmsd=unsigned(dist*precision);
	    
	    fout.write((char*)&irmsd, sizeof(unsigned));
	  }
	}
      }	
    }
  }
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
