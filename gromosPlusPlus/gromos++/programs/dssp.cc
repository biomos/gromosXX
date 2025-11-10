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
 * @file dssp.cc
 * Monitors secondary structure elements
 */

/**
 * @page programs Program Documentation
 *
 * @anchor dssp
 * @section dssp monitors secondary structure elements
 * @author ub @ref ae
 * @date 6-16-2013
 *
 * Program dssp monitors secondary structure elements for protein structures
 * over a molecular trajectory. The definitions are according to the DSSP rules
 * defined by Kabsch and Sander (Biopolymers, 22, 2577 - 2637 (1983)).  Within
 * these rules it may occur that one residue is defined as being part of two
 * different secondary-structure elements. In order to avoid duplicates in the
 * output, the following priority rules are applied: Beta Sheet/Bridge &gt;
 * 4-helix &gt; 5-helix &gt; 3-helix &gt; H-bonded turn &gt; Bend. As a
 * consequence, there may be, for instance, helices that are shorter than their
 * minimal length. 
 *
 * The program summarizes the observed occurrences of the secondary
 * structure elements and averages the different properties over the protein.
 * In addition time series for every type of secondary structure element are
 * written to file.
 *
 * NOTE: This program considers intra-molecular secondary structure only. If
 * you are interested in inter-molecular secondary protein structure, e.g. a
 * beta-sheet out of two inter-molecular beta strands, this program cannot detect
 * it.
 *
 *
 * * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gathermethod&gt;] </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms" to consider&gt; </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;]</td></tr>
 * <tr><td> [\@nthframe</td><td>&lt;write every nth frame&gt; (default is 1)] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  dssp
    @topo      ex.top
    @pbc       r
    @atoms     1:a
    @time      0 1
    @nthframe  10
    @traj      ex.tr
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/Neighbours.h"
#include "../src/utils/Hbond_calc.h"
#include "../src/utils/Dssp.h"
#include "../src/utils/groTime.h"
#include "../src/gromos/Exception.h"

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace gmath;
using namespace utils;
using namespace std;
using utils::AtomSpecifier;

int main(int argc, char **argv){
  
  Argument_List knowns;
  knowns << "topo" << "pbc" << "atoms" << "time" << "nthframe" << "traj";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo      <molecular topology file>\n";
  usage += "\t@pbc       <boundary type> [<gathermethod>]\n";
  usage += "\t@atoms     <atoms to consider>\n";
  usage += "\t[@time      <time and dt>]\n";
  usage += "\t[@nthframe <write every nth frame> (default is 1)]\n";
  usage += "\t@traj      <trajectory files>\n";
  
  
  try{
    Arguments args(argc, argv, knowns, usage);
    
    args.check("topo",1);
    InTopology it(args["topo"]);
    System sys(it.system());
    
    int nthFrame = args.getValue<int>("nthframe", false, 1);
    
    // get the protein atoms
    AtomSpecifier prot(sys);
    {
      Arguments::const_iterator iter=args.lower_bound("atoms");
      Arguments::const_iterator to=args.upper_bound("atoms");
      if (iter == to)
	throw Arguments::Exception("argument @atoms is required");
      cout << "# protein atoms considered ";
      for(;iter!=to;iter++) {
	prot.addSpecifier(iter->second.c_str());
	cout << iter->second.c_str() << " ";
      }
    }
    cout << endl;

    Dssp SecStr(sys,args);
    
    //get time
    Time time(args);
    
    SecStr.calcnumres(prot, sys);
    SecStr.determineAtoms(prot);
    SecStr.calcHintra_init(prot); 
    
    InG96 ic(0,nthFrame);
    
    for(Arguments::const_iterator 
	  iter=args.lower_bound("traj"),
	  to=args.upper_bound("traj");
	iter!=to; ++iter) {
      ic.open((iter->second).c_str());
      while(!ic.eof()) {
	ic >> sys >> time;	
      if (ic.stride_eof()) break;   
	  SecStr.calcHb_Kabsch_Sander();
	  SecStr.calc_Helices();
	  SecStr.calc_Betas();
	  SecStr.calc_Bends();
	  SecStr.filter_SecStruct();
	  SecStr.writeToFiles(time.time());
	  SecStr.keepStatistics();
      }
      ic.close();
    }
    SecStr.writeSummary(cout);
    
  }
  catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}




