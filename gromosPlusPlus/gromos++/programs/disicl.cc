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
 * @file disicl.cc
 * Dihedral-based secondary structure classification for proteins and nucleic acids
 */

/**
 * @page programs Program Documentation
 *
 * @anchor disicl
 * @section disicl secondary structure classification
 * @author @ref mp
 * @date 28-07-2015
 *
 * Program disicl classifies secondary structure elements in proteins and 
 * nucleic acids based on dihedral angles (Nagy & Oostenbrink 2014, JCIM 54(1),278 
 * and Nagy & Oostenbrink 2014, JCIM 54(1),266). 
 * 
 * Angle, region and class definitions are read from a user-specified @ref DisiclLibrary 
 * "library file". The program will warn about overlapping 
 * regions and region limits that are outside the chosen periodic range and will
 * abort if two classes have the same definition. 
 *
 * The program writes out classification statistics per residue and 
 * averaged over all residues (stat_disicl.out). Timeseries are written
 * for each class (class_XXX.dat) and for the dihedral angles (ts_disicl.dat).
 *
 * The program provides an option (pdbstride) to write out pdb files containing
 * class information in the b-factor column for visualization using the "color
 * by b-factor"-function in your favorite visualization software. An additional 
 * pdb is created (colorlegend.pdb), which can be used as a kind of legend for 
 * the class color code when loaded into the visualization software together 
 * with the output pdbs. 
 * 
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gathermethod&gt;] </td></tr>
 * <tr><td> [\@atoms</td><td>&lt;@ref utils::AtomSpecifier "atoms to include"&gt; (default: all atoms)] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> [\@stride</td><td>&lt;take every n-th frame&gt;] </td></tr>
 * <tr><td> [\@skip</td><td>&lt;skip n first frames&gt;] </td></tr>
 * <tr><td> [\@periodic</td><td>&lt;dihedral angle periodic range (default: -180 180)&gt;] </td></tr>
 * <tr><td> [\@nots</td><td>&lt;do not write timeseries of the dihedrals&gt;] </td></tr>
 * <tr><td> [\@pdbstride</td><td>&lt;write out pdb coordinates every n steps (has to be >= stride)&gt;] </td></tr>
 * <tr><td> \@lib</td><td>&lt;library file&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 * 
 *
 * Example:
 * @verbatim
  disicl
    @topo ex.top
    @pbc  r cog
    @atoms 1:a
    @time 0 2
    @lib  disicl_prot.lib
    @stride 1
    @skip 10
    @pdbstride 20
    @periodic -180 180
    @nots
    @traj ex.tr
 @endverbatim
 *
 * <hr>
 */
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <time.h>
#include <vector>
#include <iomanip>
#include <iostream>
#include <cassert>
#include <sstream>
#include <string>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/Ginstream.h"
#include "../src/gcore/System.h"
#include "../src/gio/OutCoordinates.h"
#include "../src/gio/OutPdb.h"
#include "../src/fit/PositionUtils.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/groTime.h"
#include "../src/gromos/Exception.h"
#include "../src/utils/Disicl.h"


using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;



int main(int argc, char **argv) {
  clock_t t;
  t=clock();

  Argument_List knowns;
  knowns << "topo" << "pbc" << "atoms" << "time" << "stride" << 
            "skip" << "nots" << "periodic" <<"pdbstride" << "lib" <<"traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo        <molecular topology file>\n";
  usage += "\t@pbc         <boundary type> [<gathermethod>]\n";
  usage += "\t[@atoms      <atoms to include> (default: all solute atoms)]\n";
  usage += "\t[@time       <time and dt>]\n";  
  usage += "\t[@stride     <use every n-th frame>]\n";
  usage += "\t[@skip       <skip n first frames>]\n";
  usage += "\t[@nots       <do not write dihedral timeseries>]\n";
  usage += "\t[@periodic   <dihedral angle periodic range (default: -180 180)>]\n";  
  usage += "\t[@pdbstride  <write pdb coordinates every n steps>]\n";
  usage += "\t@lib         <library file defining angles, regions and classes>\n";
  usage += "\t@traj        <trajectory files>\n";
 
  try{
    Arguments args(argc, argv, knowns, usage);

    // get the @time argument
    utils::Time time(args);
    
    bool do_tser = true;
    if (args.count("nots") >= 0)
      do_tser = false;

    //  read topology
    args.check("topo", 1);
    InTopology it(args["topo"]);
    
    System sys(it.system());
    System refSys(it.system());

    // make sure there is an argument @pbc
    args.check("pbc", 1);
    
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);
    
    // get the protein atoms
    AtomSpecifier prot(sys);
    if (args.count("atoms") > 0) {
      {
        Arguments::const_iterator iter=args.lower_bound("atoms");
        Arguments::const_iterator to=args.upper_bound("atoms");
	    cerr << "# atom selection: ";
        for(;iter!=to;iter++) {
	      prot.addSpecifier(iter->second.c_str());
	      cerr << iter->second.c_str() << " ";
        }
      }
      cerr << endl;
    }
    else {
	  prot.addSpecifier("a:a");  
	  cerr << "# selected all atoms\n";
    }
    
    Dscl dscl;
    
    // read periodic range
    dscl.setPeriodic(-180,180); // default
    vector<int> p=args.getValues<int>("periodic", 2, false,dscl.getPeriodic());
    dscl.setPeriodic(p[0],p[1]);
    cerr << "# periodic range = " << p[0] << " to " << p[1] << endl;
    
    // read library
    if (args.count("lib") == 1) {
      Ginstream lib(args["lib"]);
      cerr << "# LIBRARY = " << args["lib"] << endl;
      dscl.readLibrary(lib);
    }
    else {
      throw Arguments::Exception("@lib expects exactly one library file!");
    }
    
    //dscl.print_angles();
    //dscl.print_regions();
    //dscl.print_classes();
    
    dscl.determineAtoms(prot, sys);
    dscl.getPropSpec(sys,pbc);
    if (do_tser)
      dscl.writeHeader();

    int skip = args.getValue<int>("skip", false, 0);
    int stride = args.getValue<int>("stride", false, 1);
    
    // parse outformat
    bool write_pdb=false;
    
    int pdbstride=1;
    if (args.count("pdbstride") >= 0) {
      write_pdb=true;
      if (args.count("pdbstride") == 0)
        pdbstride=1;
      else
        pdbstride=args.getValue<int>("pdbstride",false,1);
    }
    
    // define input coordinate
    InG96 ic(skip,stride);
    
    dscl.initSummary();
    dscl.initTimeseries();
    
    unsigned int  frameNum = 0;
    unsigned int skipFrame = 0;

    // loop over all trajectories
    for (Arguments::const_iterator 
	  iter=args.lower_bound("traj"),
	  to=args.upper_bound("traj");
	  iter!=to; ++iter) {
	  // open file
      ic.open((iter->second).c_str());
      
      // loop over single trajectory
      while(!ic.eof()) {
        ic.select("SOLUTE");
        if (ic.stride_eof()) break;
        ic >> sys >> time;      

        (*pbc.*gathmethod)();        
        
        dscl.calcDih();
        dscl.classifyRegions();
        dscl.classifyClasses(time.time());
        if (do_tser)
          dscl.writeDihTs(time.time());
        dscl.keepStatistics();
        
        if (write_pdb && !skipFrame) {
          OutCoordinates *oc;
          ofstream os;
          oc = new OutPdb();
          ostringstream pdbName;
          pdbName << "out_"<< setw(5)<< setfill('0') << frameNum << ".pdb";
          os.open(pdbName.str().c_str());
          oc->open(os);
          oc->writeTitle(dscl.pdbTitle());
          dscl.getBfactorValues(sys);
          oc->writeTimestep(time.steps(), time.time());
          *oc << sys; //prot;
          os.close();
          delete oc;
          oc=0;
        }
          
        frameNum++;
        skipFrame++;
        skipFrame %= pdbstride;
        
      }
      ic.close();
    } 
    //if (write_pdb) os.close();
    dscl.writeStatistics( frameNum, do_tser);
    dscl.closeTimeseries();
    if (write_pdb) dscl.writePdbColorLegend();
  } 
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  
  t = clock() - t;
  printf ("# Runtime: %f seconds\n",((float)t)/CLOCKS_PER_SEC);
  return 0;
}

