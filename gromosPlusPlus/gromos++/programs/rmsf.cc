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
 * @file rmsf.cc
 * calculates atom-positional root-mean-square fluctuations
 */

/**
 * @page programs Program Documentation
 *
 * @anchor rmsf
 * @section rmsf atom-positional root-mean-square fluctuations
 * @author @ref mk 
 * @date 26. 7. 2006
 *
 * Program rmsf calculates atom-positional root-mean-square fluctuations 
 * (rmsf) around average positions for selected atoms over a trajectory. If requested, a
 * rotational fit to a reference structure is performed for every structure in
 * the trajectory. Different sets of atoms can be specified for the fitting 
 * procedure and for the calculation of the rmsf. If no fit is required "no"
 * should be given.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gathermethod&gt;] </td></tr>
 * <tr><td> \@atomsrmsf</td><td>&lt;@ref AtomSpecifier "atoms" to consider for rmsf&gt; </td></tr>
 * <tr><td> [\@atomsfit</td><td>&lt;@ref AtomSpecifier "atoms" to consider for fit&gt;] </td></tr>
 * <tr><td> [\@ref</td><td>&lt;reference coordinates(if absent, the first frame of \@traj is reference)&gt;] </td></tr>
 * <tr><td> [\@outpdb</td><td>&lt;write average structure in pdb format with rmsf values in the b-factor column&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
  rsmf
    @topo       ex.top
    @pbc        r
    @atomsrmsf  1:CA
    @atomsfit   1:CA,C,N
    @ref        exref.coo
    @outpdb     aver.pdb
    @traj       ex.tr
@endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <string>
#include <fstream>
#include <vector>
#include <iomanip>
#include <iostream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Matrix.h"
#include "../src/gmath/Vec.h"
#include "../src/gio/OutG96S.h"
#include "../src/gio/OutPdb.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gmath;
using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;

int main(int argc, char **argv){

  Argument_List knowns; 
  knowns <<"topo" << "traj" << "atomsfit" << "atomsrmsf" << "pbc" << "ref" << "outpdb";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo        <molecular topology file>\n";
  usage += "\t@pbc         <boundary type> [<gathermethod>]\n";
  usage += "\t@atomsrmsf   <atoms to consider for rmsf>\n";
  usage += "\t[@atomsfit   <atoms to consider for fit>]\n";
  usage += "\t[@ref         <reference coordinates(if absent, the first frame of @traj is reference)>]\n";
  usage += "\t[@outpdb      <write average structure in pdb format with rmsf values in the b-factor column>]\n";
  usage += "\t@traj        <trajectory files>\n";


  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);


  try{
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System refSys(it.system());

    // System for calculation
    System sys(refSys);
    
    // read reference coordinates...
    InG96 ic;

    try{
      args.check("ref",1);
      ic.open(args["ref"]);
    }
    catch(const Arguments::Exception &){
      args.check("traj",1);
      ic.open(args["traj"]);
    }
    ic.select("ALL");
    ic >> refSys;
    ic.close();


    // Parse boundary conditions
    //Boundary *pbc = BoundaryParser::boundary(refSys, args);
    // parse gather method
    //Boundary::MemPtr gathmethod = args::GatherParser::parse(args);


    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(refSys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

    // gather reference system
    (*pbc.*gathmethod)();
    delete pbc;
    
    // System for calculation
    //System sys(refSys);
    //get the atoms for the rmsf
    AtomSpecifier rmsfatoms(sys); 
    {
      Arguments::const_iterator iter = args.lower_bound("atomsrmsf");
      Arguments::const_iterator to = args.upper_bound("atomsrmsf");
      
      for(;iter!=to;iter++){
	string spec=iter->second.c_str();
	rmsfatoms.addSpecifier(spec);
      }
    }
    
    std::string outpdb;
    bool writepdb=false;
    if (args.count("outpdb") >= 0) {
        writepdb=true;
        outpdb="aver.pdb";
        if (args.count("outpdb") > 0) outpdb=args.getValue<std::string > ("outpdb", false, "aver.pdb");
    }

    //get the atoms for the fit
    AtomSpecifier fitatoms(refSys);
    if (args.count("atomsfit") > 0) {
      Arguments::const_iterator iter = args.lower_bound("atomsfit");
      Arguments::const_iterator to = args.upper_bound("atomsfit");

      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        fitatoms.addSpecifier(spec);
      }
    } else {
      cout << "# @atomsrmsf atoms are taken for fit." << endl;
      const vector<string> & spec = rmsfatoms.toString();

      for (vector<string>::const_iterator it = spec.begin(), to = spec.end();
              it != to; ++it) {
        fitatoms.addSpecifier(*it);
      }
    }
    
    
    // Parse boundary conditions for sys
    pbc = BoundaryParser::boundary(sys, args);
    RotationalFit * rf = NULL;
    if (fitatoms.size()) {
      Reference * ref = new Reference(&refSys);
      ref->addAtomSpecifier(fitatoms);
      rf = new RotationalFit(ref);
    }
    
    
    int numFrames = 0;
    
    //vectors to store the positions, average position and eventually the rmsf value
    
    vector<Vec> apos, firstpos;
    vector<double> apos2;
    vector<double> rmsf;
    
    //loop over trajectory
    for(Arguments::const_iterator iter=args.lower_bound("traj");
	iter!=args.upper_bound("traj"); ++iter){
      ic.open(iter->second);
      // loop over all frames
      while(!ic.eof()){
        // read frame
        ic.select("ALL");
	ic >> sys;
        // initalize after it is read
        if (numFrames == 0) {
          if (rmsfatoms.size() == 0)
            throw gromos::Exception("rmsf",
                  "No atoms specified for RMSF calculation");
          rmsfatoms.sort();
          apos.resize(rmsfatoms.size(), Vec(0.0, 0.0, 0.0));
          apos2.resize(rmsfatoms.size(), 0.0);
          rmsf.resize(rmsfatoms.size(), 0.0);
          for(unsigned int i = 0; i < rmsfatoms.size(); ++i) {
            firstpos.push_back(rmsfatoms.pos(i));
          }
        }
	numFrames++;

	
	(*pbc.*gathmethod)();
        if (fitatoms.size())
	  rf->fit(&sys);
	
	// calculate <r> and <r^2>
	for(unsigned int i=0; i< rmsfatoms.size(); ++i){
          const Vec & gathpos = rmsfatoms.pos(i); // pbc->nearestImage(firstpos[i], rmsfatoms.pos(i), sys.box());
	  apos[i] += gathpos;
	  apos2[i] += gathpos.abs2();
	}

      }
      ic.close();
    } //end loop over trajectory
    
    // calculate the rmsf's
    for(unsigned int i=0; i < rmsfatoms.size(); ++i){
      apos2[i]/=numFrames;
      apos[i]/=numFrames;
      
      rmsf[i] = sqrt(apos2[i] - apos[i].abs2());
    }
    
    //spit out results
    cout << "#\n#  at          rmsf name\n";
     
    for (unsigned int i=0; i < rmsfatoms.size(); ++i) {
      cout.precision(8);
      cout << setw(5) << i+1
	   << setw(14) << rmsf[i]
	   << ' ' << setw(4) << rmsfatoms.name(i)
	   << endl;
	   if (writepdb) sys.mol(rmsfatoms.mol(i)).setBfac(rmsfatoms.atom(i),rmsf[i]);
    }
   
    if (rf != NULL) {
      delete rf->getReference();
      delete rf;
    }
    
    if (writepdb) {
    // get a system
      for(unsigned int i=0; i< rmsfatoms.size(); ++i){
        rmsfatoms.pos(i) = apos[i];
      }
      ofstream fout(outpdb.c_str());
      OutPdb oc(fout);
      oc.writeTitle("Average structure\n");
      oc << rmsfatoms;
    }
    
  }
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


