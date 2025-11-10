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
 * @file mdf.cc
 * minimum distance function
 */

/**
 * @page programs Program Documentation
 *
 * @anchor mdf
 * @section mdf minimum distance function
 * @author @ref co
 * @date 31.7.2006
 *
 * Program mdf calculates and lists, for a given, central set of atoms, the distance to the
 * nearest atom belonging to a second set of atoms. For every selected atom, an
 * output file is written with the minimum distance to, and an atom specifier
 * for, the nearest atom. This program also works for @ref VirtualAtom "virtual atoms".
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;]</td></tr>
 * <tr><td>[\@time</td><td>&lt;@ref utils::Time "time dt"&gt;]</td></tr>
 * <tr><td> \@centre</td><td>&lt;central group of @ref AtomSpecifier "atoms"&gt; </td></tr>
 * <tr><td> \@with</td><td>&lt;group of @ref AtomSpecifier "atoms" from which to find the nearest atom&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
  mdf
    @topo     ex.top
    @pbc      r
    @time     0 0.5
    @centre   1:3-5
    @with     s:OW
    @traj ex.tr
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cfloat>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/bound/Boundary.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Solvent.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
#include "../src/gromos/Exception.h"

using namespace gcore;
using namespace gmath;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;
using namespace std;


int main(int argc, char **argv){

  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "centre" << "with" << "nsm" << "traj";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@pbc    <boundary type> [<gather method>]\n";
  usage += "\t@time   <time and dt>\n";
  usage += "\t@centre <atoms to take as centre>\n";
  usage += "\t@with   <atoms to calculate the distance for>\n";
  usage += "\t@traj   <trajectory files>\n";
  
 
try{
  Arguments args(argc, argv, knowns, usage);

  // read topology
  args.check("topo",1);
  InTopology it(args["topo"]);
  System sys(it.system());

  //   get simulation time
  Time time(args);

  // set centre atoms
  AtomSpecifier centre(sys);
  
  {
    Arguments::const_iterator iter=args.lower_bound("centre");
    Arguments::const_iterator to=args.upper_bound("centre");
    for(; iter!=to; ++iter){
      centre.addSpecifier(iter->second);
    }
  }
  // set atom to consider
  AtomSpecifier with(sys);
  {
    Arguments::const_iterator iter=args.lower_bound("with");
    Arguments::const_iterator to=args.upper_bound("with");
    for(; iter!=to; ++iter){
      with.addSpecifier(iter->second);
    }
  }
  
  // Parse boundary conditions
  Boundary *pbc = BoundaryParser::boundary(sys, args);

  // define input coordinate
  InG96 ic;

  // loop over all trajectories
  int count_frame=0;
  
  // open centre.size() files
  vector<ofstream *> fout(centre.size());
  
  for(unsigned int i=0; i<centre.size(); i++){
    ostringstream os;
    if(centre.mol(i)!=-3)
      os << "MIN_" << centre.mol(i)+1 << ":" << centre.atom(i)+1 
	 << ".dat";
    else
      os << "MIN_va_" << i+1 <<".dat";
    
    fout[i] = new ofstream(os.str().c_str());
    
  }
  
  for(Arguments::const_iterator 
        iter=args.lower_bound("traj"),
        to=args.upper_bound("traj");
      iter!=to; ++iter){

    // open file
    ic.open((iter->second).c_str());
    ic.select("ALL");
    
   
    // loop over single trajectory
    while(!ic.eof()){
      ic >> sys >> time;
      pbc->gather();
      
      // loop over the centre atoms
      int start=0;
      int minat=0;
      
      // now really loop over the centre atoms
      for(unsigned int i=start;i<centre.size();i++){
//	double min2=sys.box().K().abs()*sys.box().K().abs();
//      The box dimension might not be the appropriate initial distance.
//      Therefore simply use the maximum finite representable floating-point number:
        double min2=DBL_MAX;
	//loop over the atoms to consider
        for(unsigned int j=0;j<with.size();j++){
          //calculate distance only if this atom is not the current centre
          if(!(with.mol(j)==centre.mol(i)&&with.atom(j)==centre.atom(i))){
	    Vec tmp=pbc->nearestImage(centre.pos(i),
				      with.pos(j),
				      sys.box());
	    tmp-=centre.pos(i);
	    
            if(tmp.abs2()<min2) {
              min2=tmp.abs2();
              minat=j;
	    }
	    
	  }
	}
	//write out min dist
        (*fout[i]) << time << "\t" << sqrt(min2) << "\t# " << with.toString(minat) << endl;
	
      }
      count_frame++;
    }
    ic.close();
  }
  
  //close the files
  for(unsigned int i=0;i<centre.size();i++){
    fout[i]->close();
    delete fout[i];
    fout[i] = NULL;
  }
  
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}




























