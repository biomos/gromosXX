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
 * @file follow.cc
 * Display a trajectory of selected atoms in 3D
 */

/**
 * @page programs Program Documentation
 *
 * @anchor follow
 * @section follow Display a trajectory of selected atoms in 3D
 * @author @ref co
 * @date 31-10-08
 *
 * Program follow can create a 3D trace of selected atoms through time. The 
 * program always takes the nearest image with respect to the previous position
 * of the particle. For every atom that is selected, a pdb file is written out
 * (FOLLOW_x.pdb) in which the trajectory is indicated in the CONECT entries. 
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;topology&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> \@dim</td><td>&lt;dimensions to consider&gt; </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms" to follow&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  follow
    @topo  ex.top
    @pbc   t
    [@time  0 1]
    @dim   x y z
    @atoms 1:5
    @traj  ex.tr
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>


#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Solvent.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
#include "../src/gromos/Exception.h"

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

using namespace std;

string molecule(int m, int a);
string dimension(int dim);
void writepdb(ofstream &out, int count, int mol, string atname, 
	      Vec v, int dim[3], int ndim);
void writeCON(ofstream &out, int count, int jump, int start);


int main(int argc, char **argv){

  Argument_List knowns; 
  knowns << "topo" << "pbc" << "time" << "dim" << "atoms" << "traj";

  string usage = argv[0];
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t[@time   <time and dt>]\n";
  usage += "\t@dim    <dimensions to consider>\n";
  usage += "\t@atoms  <atoms to follow>\n";
  usage += "\t@traj   <trajectory files>\n";
  
 
try{
  Arguments args(argc, argv, knowns, usage);

  // get the @time argument
    utils::Time time(args);

  // get the relevant dimensions
  int ndim=3;
  int dim[3]={0,1,2};
  {
    Arguments::const_iterator iter=args.lower_bound("dim");
    if(iter!=args.upper_bound("dim")){
      ndim=0;
      
      string dum=iter->second.c_str();
      if(dum=="x")      { dim[ndim]=0; ndim++;}
      else if(dum=="y") { dim[ndim]=1; ndim++;}
      else if(dum=="z") { dim[ndim]=2; ndim++;}
      iter++;
    }
    if(iter!=args.upper_bound("dim")){
      string dum=iter->second.c_str();
      if(dum=="x")      { dim[ndim]=0; ndim++;}
      else if(dum=="y") { dim[ndim]=1; ndim++;}
      else if(dum=="z") { dim[ndim]=2; ndim++;}
      iter++;
    }
    if(iter!=args.upper_bound("dim")){
      string dum=iter->second.c_str();
      if(dum=="x")      { dim[ndim]=0; ndim++;}
      else if(dum=="y") { dim[ndim]=1; ndim++;}
      else if(dum=="z") { dim[ndim]=2; ndim++;}
    }
  }
  
  //  read topology
  args.check("topo",1);
  InTopology it(args["topo"]);
  System sys(it.system());
  
  // parse boundary conditions
  Boundary *pbc = BoundaryParser::boundary(sys, args);

  // set atom number
  AtomSpecifier at(sys);
  {
    Arguments::const_iterator iter=args.lower_bound("atoms");
    Arguments::const_iterator to=args.upper_bound("atoms");
    for(;iter!=to;iter++)
      at.addSpecifier(iter->second.c_str());
  }
  // we need to store the old coordinates of the atoms to follow
  vector<Vec> oldpos;
  
  for(unsigned int i=0; i<at.size(); i++){
    oldpos.push_back(Vec(0.0,0.0,0.0));
  }
  // open at.size() files to write the trajectories to pdb
  vector<ofstream *> opdb(at.size());

  for(unsigned int i=0; i<at.size(); i++){
    stringstream os;
    string s=molecule(at.mol(i), at.atom(i));
    s=s.substr(0,s.find(' '));
    
    os << "FOLLOW_" << s << ".pdb" << ends;

    opdb[i] = new ofstream(os.str().c_str());
    (*opdb[i]) << "TITLE following coordinates of " 
	       << molecule(at.mol(i), at.atom(i)) << endl;
  }
  
    
  // print title
  cout << setw(6) << "#     ";
  for(unsigned int i=0; i<at.size(); i++)
    cout << setw(ndim*12) << molecule(at.mol(i), at.atom(i));
  cout << endl;
  cout << setw(6) << "# time";
  for(unsigned int i=0; i<at.size(); i++)
    for(int k=0; k<ndim; k++)
      cout << setw(12) << dimension(dim[k]);
  cout << endl;
  
  
     
  int frames=1;
  int count=1;
  
  InG96 ic;
  
  // loop over all trajectories
  for(Arguments::const_iterator 
      iter=args.lower_bound("traj"), to=args.upper_bound("traj");
      iter!=to; ++iter){

      // open file
    ic.open((iter->second).c_str());
    ic.select("ALL");
    
      // loop over single trajectory
    while(!ic.eof()){
      ic >> sys >> time;

     // if(count==1){
	//for(int i=0; i<at.size(); i++){
	//  oldpos[i]=at.pos(i);
	//}
      //}
      
      cout << time;
      
      // loop over the atoms to consider
      for(unsigned int i=0; i<at.size(); i++){
	// gather with respect to its old position
	*at.coord(i)=
	  pbc->nearestImage(oldpos[i], *at.coord(i), sys.box());
	// print out the relevant coordinates
	for(int j=0; j<ndim; j++){
	  cout << setw(12) << (*at.coord(i))[dim[j]];
	}
	writepdb(*opdb[i], count, i, at.name(i), *at.coord(i), dim, ndim);
	
	// copy the current system to oldsys
	oldpos[i]=*at.coord(i);
        count++;
	
      }
      cout << endl;
      
      frames++;
    }
    
    ic.close();
  }
  for(unsigned int i=0; i<at.size(); i++){
    
    writeCON(*opdb[i], count, at.size(), i);

    opdb[i]->close();
    delete opdb[i];
    opdb[i] = NULL;
    
  }
  
  
}
  
 
  
catch (const gromos::Exception &e){
  cerr << e.what() << endl;
  exit(1);
}
return 0;
}

string molecule(int m, int a)
{
  stringstream os;
  if(m<0) os << "s";
  else os << m+1;
  os << ":" << a+1;
  if(m<0) os << " (solv)";
  else os << " (solu)";
  os << ends;
  return os.str().c_str();
}


string dimension(int dim)
{
  if(dim==0) return "x";
  if(dim==1) return "y";
  if(dim==2) return "z";
  return "?";
  
}

void writepdb(ofstream &out, int count, int mol, string atname, 
	      Vec v, int dim[3], int ndim)
{
  out.setf(ios::fixed, ios::floatfield);
  out.setf(ios::unitbuf);
  out.precision(3);

  out << "HETATM";
  out.setf(ios::right, ios::adjustfield);
  out << setw(5) << count;
  
  out.setf(ios::left, ios::adjustfield);
  out << "  " <<setw(4) << atname.c_str();
  out << setw(4) << "FLW  " << setw(4) << mol+1 << "    ";
  out.setf(ios::right, ios::adjustfield);
  Vec v2(0.0,0.0,0.0);
  for(int i=0; i<ndim; i++)
    v2[dim[i]]=v[dim[i]];
  
  out  << setw(8) << v2[0]*10
       << setw(8) << v2[1]*10
       << setw(8) << v2[2]*10
       << "  1.00  0.00" << endl;
  
}

void writeCON(ofstream &out, int count, int jump, int start)
{
  for(int i=1; i<count-jump; i+=jump)
    out << "CONECT"
	<< setw(5) << i+start << setw(5) << i+jump+start << endl;
  out << "TER" << endl;
}


