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
 * @file dipole.cc
 * Calculates the dipole moment for specified atoms
 */

/**
 * @page programs Program Documentation
 *
 * @anchor dipole
 * @section dipole Calculates the dipole moment for specified atoms
 * @author @ref mk
 * @date 11-6-07
 *
 * Program dipole will calculate and print the dipole moment of molecules. By
 * default, the program will take all solute atoms into account, but the user 
 * can also specify a separate set of atoms. The dipole moment of a set of 
 * atoms carrying a net-charge is ill-defined and depends on the position of 
 * the origin. For these cases, the program allows the user to move the 
 * centre-of-geometry of the atoms to the origin.
 *
 * <B>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> [\@atoms</td><td>&lt;@ref AtomSpecifier "atoms" to include&gt; (default all solute)] </td></tr>
 * <tr><td> [\@cog</td><td>(move molecule to centre-of-geometry)] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  dipole
    @topo  ex.top
    @pbc   r
    [@time  0 0.2]
    @atoms 1:1-30
    @cog
    @traj  ex.tr
 @endverbatim
 *
 * <hr>
 */

#include <cstdlib>
#include <string>
#include <iomanip>
#include <iostream>
#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace fit;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;

int main(int argc, char **argv){

  Argument_List knowns; 
  knowns << "topo" << "pbc" << "time" << "atoms" << "cog" << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo    <molecular topology file>\n";
  usage += "\t@pbc     <boundary type> [<gather method>]\n";
  usage += "\t[@time    <time and dt>]\n";
  usage += "\t[@atoms   <atoms to include> (default all solute)]\n";
  usage += "\t[@cog    (move molecule to centre-of-geometry)]\n";
  usage += "\t@traj    <trajectory files>\n";
  
 
  try{
    Arguments args(argc, argv, knowns, usage);

    // get the @time argument
    utils::Time time(args);

    //  read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    System refSys(it.system());

    //get the atoms to be included
    utils::AtomSpecifier atoms(sys);
    if(args.count("atoms")>0)
      for(Arguments::const_iterator iter=args.lower_bound("atoms"), 
            to=args.upper_bound("atoms"); iter!=to; ++iter){
        atoms.addSpecifier(iter->second);
      }
    else
      for(int m=0; m<sys.numMolecules(); m++)
        for(int a=0; a<sys.mol(0).numAtoms(); a++)
          atoms.addAtom(m,a);
 
 
    //determine net charge
    double ncharge=0;
    double nchargepa=0;
    for(unsigned int i=0; i < atoms.size(); i++){
      ncharge+=atoms.charge(i);
    }
    nchargepa=ncharge/atoms.size();
    
    // Move to cog?
    bool lcog=false;
    if(args.count("cog") >= 0) lcog=true;
    
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);
    
    // write a title
    cout << "#\n";
    if(ncharge!=0.0)
      cout << "# WARNING the specified atoms carry a net charge ( "
	   << ncharge << ")\n"
	   << "#         this means that the dipole depends on the origin\n"
	   << "#\n";
    cout << "# molecule is ";
    if(!lcog) cout << "not ";
    cout << "translated to its centre-of-geometry before dipole calculation\n";
    if(lcog && ncharge==0.0) 
      cout << "# (this is irrelevant as the netcharge is zero)\n";
    cout << "#\n";
    cout << "#     time  magnitude          x          y          z\n";

    // define input coordinate
    InG96 ic;
    
    // variables for averaging
    Vec tot_dip(0.0,0.0,0.0);
    double tot_mag=0.0;
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
	ic >> sys >> time;
	(*pbc.*gathmethod)();
	Vec cog(0.0,0.0,0.0);
	if(lcog){
	  for(unsigned int i=0; i<atoms.size(); ++i)
	    cog+=atoms.pos(i);
	  cog/=atoms.size();
	}
	Vec dipole(0.0,0.0,0.0);
	for(unsigned int i=0; i<atoms.size(); ++i){
	  dipole += (atoms.charge(i)-nchargepa)*(atoms.pos(i)-cog);
	}
	
	cout << time << " ";
	cout << setw(10) << dipole.abs() << " " 
	     << setw(10) << dipole[0] << " "
	     << setw(10) << dipole[1] << " "
	     << setw(10) << dipole[2] <<  endl;
	tot_mag+=dipole.abs();
	tot_dip+=dipole;
	numFrames++;
      }
      ic.close();
    }
    cout << "#\n"
	 << "# averages "
	 << setw(10) << tot_mag/numFrames << " "
	 << setw(10) << tot_dip[0]/numFrames << " "
	 << setw(10) << tot_dip[1]/numFrames << " "
	 << setw(10) << tot_dip[2]/numFrames
	 << endl;
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
