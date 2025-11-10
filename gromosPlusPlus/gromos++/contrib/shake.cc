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
 * @file shake.cc
 * applies the SHAKE algorithm to a trajectory.
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor shake
 * @section shake SHAKE a trajectory
 * @author @ref mk @ref ns
 * @date 19. 2. 2008
 *
 * Program shakes applies the SHAKE algorithm to solute structures in a 
 * trajectory. 
 *
 * This program needs the GROMOS XX library to run.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gathermethod&gt;] </td></tr>
 * <tr><td> [\@tol</td><td>&lt;SHAKE tolerance: default 0.0001&gt;] </td></tr>
 * <tr><td> [\@atomsfit</td><td>&lt;@ref AtomSpecifier "atoms" to consider for fit&gt;] </td></tr>
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
    @traj       ex.tr
@endverbatim
 *
 * <hr>
 */

#include "../config.h"

#ifdef HAVE_MDPP

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <sstream>
#include <list>

// the md++ stuff

#include <md++/math/gmath.h>
#include <md++/util/debug.h>
#include <md++/util/error.h>
#include <md++/util/timing.h>

#include <md++/io/argument.h>
#include <md++/io/message.h>


#include <md++/algorithm/algorithm.h>
#include <md++/topology/topology.h>
#include <md++/simulation/simulation.h>
#include <md++/math/fft.h>
#include <md++/configuration/configuration.h>

#include <md++/algorithm/algorithm/algorithm_sequence.h>

#include <md++/interaction/interaction.h>
#include <md++/interaction/interaction_types.h>
#include <md++/interaction/forcefield/forcefield.h>

#include <md++/util/parse_verbosity.h>

#include <md++/io/topology/in_topology.h>

#include <md++/algorithm/constraints/create_constraints.h>
#include <md++/util/create_simulation.h>

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

int main(int argc, char **argv){

  Argument_List knowns;
  knowns << "topo" << "pbc" << "tol" << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t[@tol   <tolerance: default 0.0001>]\n";  
  usage += "\t@traj   <trajectory files>\n";
  
 
  try{
    Arguments args(argc, argv, knowns, usage);
    
    //   get SHAKE tolerance
    double tol = 0.0001; 
    {
      Arguments::const_iterator iter=args.lower_bound("tol");
      if(iter!=args.upper_bound("tol"))
	tol=atof(iter->second.c_str());
    }
      
    //  read topology
    args.check("topo",1);
    InTopology it(args["topo"]);
    System sys(it.system());

    //==================================================
    //=== XX INITIALIZATION                           ==
    //==================================================
    util::simulation_struct a_xx_sim;

    a_xx_sim.sim.param().system.npm = 1;
    a_xx_sim.sim.param().system.nsm = 0;
    a_xx_sim.sim.param().constraint.ntc = 3;
    a_xx_sim.sim.param().constraint.solute.algorithm = simulation::constr_shake;
    a_xx_sim.sim.param().constraint.solvent.algorithm = simulation::constr_off;
    a_xx_sim.sim.param().constraint.solute.shake_tolerance = tol;
    a_xx_sim.sim.param().constraint.solvent.shake_tolerance = tol;
    
    {
      // create a XX In_Topology
      io::In_Topology in_topo;
      
      // no output...
      in_topo.quiet = true;
    
      if (util::create_simulation(args["topo"],
				  "",
				  "",
				  "",
				  a_xx_sim,
				  in_topo)){
	std::cerr << "creating the XX system failed!" << std::endl;
	return 1;
      }
    
      if (algorithm::create_constraints(a_xx_sim.md,
					a_xx_sim.topo,
					a_xx_sim.sim,
					in_topo,
					true)){
	std::cerr << "creating the constraints algorithm failed!" << std::endl;
	return 1;
      }
      
    } // don't need the In_Topology any more...
    
    a_xx_sim.conf.resize(a_xx_sim.topo.num_atoms());

    if (args["pbc"] == "r")
      a_xx_sim.conf.boundary_type = math::rectangular;
    else if(args["pbc"] == "c")
      a_xx_sim.conf.boundary_type = math::triclinic;
    else if (args["pbc"] == "t")
      a_xx_sim.conf.boundary_type = math::truncoct;
    else if (args["pbc"] == "v")
      a_xx_sim.conf.boundary_type = math::vacuum;
    else
      throw gromos::Exception("pbc", "wrong boundary condition (only r and v)");

    // initialise arrays but do not gather
    a_xx_sim.conf.init(a_xx_sim.topo,
                       a_xx_sim.sim.param(),
                       false);
      
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    
    // set the boundary type also in md++
    switch(pbc->type()) {
      case 'r' : 
        a_xx_sim.conf.boundary_type = math::rectangular;
        break;
      case 't' :
        a_xx_sim.conf.boundary_type = math::truncoct;
        break;
      case 'c' :
        a_xx_sim.conf.boundary_type = math::triclinic; 
        break;
      case 'v' :
      default :
        a_xx_sim.conf.boundary_type = math::vacuum;
        break;
    }    
    // parse gather method
    System refSys(it.system());
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);
      
    // define input coordinate
    InG96 ic;
    OutG96 oc(cout);
    
    // get rid of messages
    io::messages.clear();
    
    // loop over all trajectories
    unsigned int num_frames = 0;
    for(Arguments::const_iterator 
	  iter=args.lower_bound("traj"),
	  to=args.upper_bound("traj");
	iter!=to; ++iter){

      // open file
      ic.open((iter->second).c_str());

      // title
      std::ostringstream is;
      is << ic.title() << endl << "shaken by GromosXX"
         << "using topology: " << args["topo"] << endl;
      oc.writeTitle(is.str());
      
      // loop over single trajectory
      while(!ic.eof()){
	ic >> sys;
	(*pbc.*gathmethod)();
        num_frames++;

	// parse the solute positions over
	for(int m=0, t=0; m < sys.numMolecules(); ++m){
	  for(int a=0; a < sys.mol(m).numAtoms(); ++a, ++t){
	    for(int d=0; d<3; ++d){
	      a_xx_sim.conf.current().pos(t)(d) = sys.mol(m).pos(a)[d];
	    }
	  }
	}
	// exchange them (for shake, we need new and old)
	a_xx_sim.conf.exchange_state();
	a_xx_sim.conf.current().pos = a_xx_sim.conf.old().pos;

	// parse the box over
        switch(a_xx_sim.conf.boundary_type) {
          case math::truncoct :
          case math::rectangular :
            a_xx_sim.conf.current().box(0) = math::Vec(sys.box().K()[0], 0.0, 0.0);
            a_xx_sim.conf.current().box(1) = math::Vec(0.0, sys.box().L()[1], 0.0);
            a_xx_sim.conf.current().box(2) = math::Vec(0.0, 0.0, sys.box().M()[2]);
            break;
          case math::triclinic :
            a_xx_sim.conf.current().box(0) = math::Vec(sys.box().K()[0],
                                                       sys.box().K()[1],
                                                       sys.box().K()[2]);
            a_xx_sim.conf.current().box(1) = math::Vec(sys.box().L()[0],
                                                       sys.box().L()[1],
                                                       sys.box().L()[2]);
            a_xx_sim.conf.current().box(2) = math::Vec(sys.box().M()[0],
                                                       sys.box().M()[1],
                                                       sys.box().M()[2]);
            break;
          case math::vacuum : // do nothing
          default : 
            break;
        }
	
	// shake
	int error = a_xx_sim.md.run(a_xx_sim.topo, a_xx_sim.conf, a_xx_sim.sim);
        if (error) {
          cerr << "Warning: SHAKE failure at frame " << num_frames << endl;
          io::messages.display();
          io::messages.clear();
        }

	// recover positions
	// parse the positions over
	for(int m=0, t=0; m < sys.numMolecules(); ++m){
	  for(int a=0; a < sys.mol(m).numAtoms(); ++a, ++t){
	    for(int d=0; d<3; ++d){
	      sys.mol(m).pos(a)[d] = a_xx_sim.conf.current().pos(t)(d);
	    }
	  }
	}
        
        // write out the system
	oc << sys;
      }	
    }

    ic.close();
    oc.close();
    
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

// no MD++
#else

#include <iostream>

int main()
{
  std::cout << "\nconfigure could not find the MD++ libraries" << std::endl
	    << "needed to run this program." << std::endl << std::endl
	    << "You need to add them to your CPPFLAGS, CXXFLAGS, LDFLAGS" << std::endl
            << "or run ./configure --with-mdpp=<path>" << std::endl << std::endl
	    << "Reconfigure and recompile to use this program" << std::endl;
  return 1;
}

#endif
