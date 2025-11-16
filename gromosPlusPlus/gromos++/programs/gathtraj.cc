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
 * @file gathtraj.cc
 * Gather a trajectory
 */

/**
 * @page programs Program Documentation
 *
 * @anchor gathtraj
 * @section gathtraj writes a gathered trajectory
 * @author @ref vk
 * @date 01.03.2002
 *
 *
 * Program gathtraj applies the periodic boundary conditions to a coordinate
 * trajectory and writes the gathered trajectory.
 *
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; &lt;gather method&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  gathtraj
    @topo        ex.top
    @pbc         r
    @traj        ex.tr
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cstdlib>
#include <string>
#include <iostream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gio/OutCoordinates.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;



int main(int argc, char **argv){

  Argument_List knowns; 
  knowns << "topo" << "pbc" << "traj";

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@traj <trajectory files>\n";
 
  try{
    Arguments args(argc, argv, knowns, usage);

    //  read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    System refSys(it.system());
    
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);

    // GatherParser
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

    
    // define input coordinate
    InG96 ic;

    // output
    OutCoordinates *oc;
    oc = new OutG96();
    oc->open(cout);  
    oc->writeTitle("gathered trajectory");

    // loop over all trajectories
    for(Arguments::const_iterator 
        iter=args.lower_bound("traj"),
        to=args.upper_bound("traj");
	iter!=to; ++iter){
      
      // open file
      ic.open((iter->second).c_str());
      
      // loop over single trajectory
      while(!ic.eof()){
	ic >> sys;
	(*pbc.*gathmethod)();
	
	*oc << sys;
      }
      
      ic.close();
      oc->close();
    }
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
