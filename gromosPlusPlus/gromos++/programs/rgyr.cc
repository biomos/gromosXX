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
 * @file rgyr.cc
 * calculates radius of gyration
 */

/**
 * @page programs Program Documentation
 *
 * @anchor rgyr
 * @section rgyr calculated radius of gyration
 * @author @ref mk 
 * @date 28. 7. 2006
 *
 * Program rgyr calculates the radius of gyration (@f$R_{gyr}@f$) for a 
 * selected set of atoms over the trajectory according to
 * @f[ R_{gyr} = \sqrt{ \frac{1}{N}\sum{(r_i - r_{com})^2}} @f]
 * where N is the number of specified atoms, @f$r_i@f$ is the position of 
 * particle i and @f$r_{com}@f$ is the centre-of-mass of the atoms. 
 * Alternatively, the radius of gyration can be calculated in a mass-weighted 
 * manner:
 * @f[ R_{gyr} = \sqrt{ \frac{1}{M}\sum{m_{i} (r_i - r_{com})^2}} @f]
 * where @f$M@f$ is the total mass of the specified atoms and @f$m_i@f$ is the 
 * mass of particle i.
 * Please note that in case atoms from more than one molecule has been chosen,
 * care should be taken in the choice of gather method to ensure a proper
 * calculation of the centre-of-mass.
 * 
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gathermethod&gt;] </td></tr>
 * <tr><td> \@time</td><td>&lt;@ref utils::Time "time and dt"&gt; </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms" to consider&gt; </td></tr>
 * <tr><td> [\@massweighted</td><td>(use massweighted formula)]</td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  rgyr
    @topo           ex.top
    @pbc            r cog
    @time           0 0.5
    @atoms          1:a
    @massweighted
    @traj           ex.tr
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include <iomanip>
#include <iostream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/groTime.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

int main(int argc, char **argv){
  
  Argument_List knowns; 
  knowns << "topo" << "pbc" << "time" << "atoms" << "massweighted" << "traj";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo           <molecular topology file>\n";
  usage += "\t@pbc            <boundary type> [<gathermethod>]\n";
  usage += "\t@time           <time and dt>\n";
  usage += "\t@atoms          < atoms to consider>\n";
  usage += "\t[@massweighted  (use massweighted formula)]\n";
  usage += "\t@traj           <trajectory files>\n";
  
  
  try{
    Arguments args(argc, argv, knowns, usage);
    
    //   get simulation time
    Time time(args);
    
    //  read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    System refSys(it.system());
    
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);
    
    // get molecules
    AtomSpecifier atom(sys);
    {
      Arguments::const_iterator iter=args.lower_bound("atoms");
      Arguments::const_iterator to=args.upper_bound("atoms");
      for(;iter!=to;iter++) {
        atom.addSpecifier(iter->second.c_str());
      }
    }
    cout << "# Radius of Gyration for: \n# ";  
    std::vector<std::string> names=atom.toString();
    for(unsigned int i=0; i< names.size(); ++i)
      cout << " " << names[i];
    cout << endl;  

    // calculate total mass
    double totalMass=0.0;
    for(unsigned int i=0; i< atom.size(); ++i)
      totalMass+=atom.mass(i);
    
    bool massweighted = false;
    if(args.count("massweighted") >=0) massweighted = true;
    
    // define input coordinate
    InG96 ic;

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
	
	//calculate cm, rgyr
	Vec cm(0.0,0.0,0.0);
	for (unsigned int i=0;i < atom.size(); i++) {
	  cm += atom.pos(i) * atom.mass(i);
	}
	cm /= totalMass;
	
	double rg=0; 

	if(massweighted) {
	  for(unsigned int i=0; i < atom.size(); ++i){
	    rg += atom.mass(i)*(atom.pos(i) - cm).abs2();
	  }
	  rg = sqrt(rg/totalMass);
	}
	else {
	  for (unsigned int i=0;i < atom.size(); i++) {
	    // should we correct for periodicity?
	    // Only if atom contains different molecules.
	    // But then cm would be wrong already. The user should just use 
	    // gathermethod cog.
	    rg += (atom.pos(i)-cm).abs2();
	  }              
	  rg = sqrt(rg/(atom.size()));
	}
	
	cout << time;
	cout << setw(15) << rg << "\n";
      }
      ic.close();
    }
  }
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
  
}

