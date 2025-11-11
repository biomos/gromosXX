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
 * @file espmap.cc
 * Calculates the vacuum electrostatic potential around a group of atoms
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor espmap
 * @section espmap Calculates the vacuum electrostatic potential around a group of atoms
 * @author @ref mk
 * @date 26-06-07
 *
 * Program espmap calculates the vacuum electrostatic potential around a user
 * specified group of atoms. It uses the atomic partial charges as defined in 
 * the topology and calculates the potential on a grid. The results are written
 * to a .pl file that can be converted to an .plt file which can be read in by
 * Gopenmol.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms" to consider&gt; </td></tr>
 * <tr><td> \@grspace</td><td>&lt;grid spacing (default: 0.2 nm) </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  espmap
    @topo     ex.top
    @pbc      r
    @atoms    a:a
    @grspace  0.2
    @traj     exref.coo
 @endverbatim
 *
 * <hr>
 */

#include <algorithm>
#include <cassert>
#include <vector>
#include <iostream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/gcore/Box.h"
#include "../src/gio/OutPdb.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/gio/OutCoordinates.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;

int main(int argc, char **argv){
  
  Argument_List knowns;
  knowns << "topo" << "pbc" << "atoms" << "grspace" << "traj";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo    <molecular topology file>\n";
  usage += "\t@pbc     <boundary type>\n";
  usage += "\t@atoms   <atoms to consider>\n";
  usage += "\t@grspace <grid spacing (default: 0.2 nm)\n";
  usage += "\t@traj    <trajectory files>\n";
  
  
  try{
    Arguments args(argc, argv, knowns, usage);
    
    //  read topology
    InTopology it(args["topo"]);
    //make system
    System sys(it.system());
    
    // which atoms considered?
    utils::AtomSpecifier atoms(sys);
    for(Arguments::const_iterator it=args.lower_bound("atoms");
	it!=args.upper_bound("atoms"); ++it){
      atoms.addSpecifier(it->second);
    }

    // get grid spacing
    double space = args.getValue("grspace", false, 0.2);
    
    System refSys(it.system());
 
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);
    
    // define input coordinate
    InG96 ic;
    Vec center (0.0, 0.0, 0.0), cog (0.0,0.0,0.0), bx(0.0,0.0,0.0),
      grmin(0.0,0.0,0.0), grmax(0.0,0.0,0.0);
    vector<Vec> espgrid;
    vector<double> esp;
    double esptmp=0;
    int nx=0,ny=0,nz=0, numFrames=0;
    double ONE_4PI_EPS0 = 138.9354;
    
    OutCoordinates *oc;
    oc = new OutPdb();
    string ext = ".pdb";
    string extpl = ".pl";
    
    
    // loop over all trajectories
    for(Arguments::const_iterator 
	  iter=args.lower_bound("traj"),
	  to=args.upper_bound("traj");
	iter!=to; ++iter){
      
      // open file
      ic.open((iter->second).c_str());
      
      // loop over single trajectory
      while(!ic.eof()){
	numFrames+=1;
	ic >> sys;
	(*pbc.*gathmethod)();
	
	// calculate the center of geometry of the aoms
	for(int i=0; i< atoms.size(); ++i)
	  cog += atoms.pos(i);
	cog /= atoms.size();
	
	//build grid positions
	Box box = sys.box();
	bx[0]= rint(box.K()[0]);bx[1]= rint(box.L()[1]);bx[2]= rint(box.M()[2]);
        //bx[0]= rint(box[0]);bx[1]= rint(box[1]);bx[2]= rint(box[2]);
	
	nx=int(rint(bx[0]/space)); ny=int (rint(bx[1]/space)); nz=int (rint(bx[2]/space));
	//make sure we have equal number of points in x,y,z
	//dont know if this is nescessary for the programs that read in grids
	//i guess it wont hurt either...
	nx=max(nx,ny); nx=max(nx,nz);
	ny=nx; nz=ny;
	
	Vec start = cog-bx/2;
	
	
	for(int i=0;i<nx;i++){
	  for(int j=0;j<ny;j++){
	    for(int k=0;k<nz;k++){ 
	      Vec tmp;      
	      tmp[0]=start[0]+space*i;
	      tmp[1]=start[1]+space*j;
	      tmp[2]=start[2]+space*k;
	      espgrid.push_back(tmp);
	      
	    }
	  }
	}
	grmin=start*10;
	grmax=espgrid[espgrid.size()-1]*10;
	
	
	//calculate esp from gromos charges
	
	for (int i=0;i < int (espgrid.size()); ++i){
	  Vec tmp = espgrid[i];
	  for(int ii=0; ii < atoms.size(); ++ii){
	    Vec rvec = atoms.pos(ii) - tmp;
	    double r = rvec.abs();
	    esptmp += ONE_4PI_EPS0*atoms.charge(ii)/r;
	  }
	  esp.push_back(esptmp);
	  esptmp=0;
	}
	
	
	char outFile[]="FRAME";
	ostringstream out;
	ostringstream outpl;
	if (numFrames < 10){
	  out << outFile <<"_"<<"0000"<< numFrames<<ext;
	  outpl << outFile <<"_"<<"0000"<< numFrames<<extpl;
	}
	else if (numFrames < 100){
	  out << outFile<<"_"<<"000"<< numFrames<<ext;
	  outpl << outFile<<"_"<<"000"<< numFrames<<extpl;
	}
	else if (numFrames < 1000){
	  out << outFile <<"_"<<"00"<< numFrames<<ext;
	  outpl << outFile <<"_"<<"00"<< numFrames<<extpl;
	}
	else {
	  out << outFile <<"_"<<"0"<< numFrames<<ext;
	  outpl << outFile <<"_"<<"0"<< numFrames<<extpl;
	}
	//write pdb
	ofstream os(out.str().c_str()); 
	oc->open(os);       
	oc->select("SOLUTE");
	oc->writeTitle(out.str());
	
	*oc << sys;
	os.close();
	
	//write pl-file
	ofstream opl;
	opl.open(outpl.str().c_str());
	
	opl << "3" << ' ' << "3" << endl;
	opl << nx << ' ' << ny << ' ' << nz << endl;
	//               opl.setf(scientific, floatfield);
	opl << grmin[2] << ' ' << grmax[2] << ' ' 
	    << grmin[1] << ' ' << grmax[1] << ' '
	    << grmin[0] << ' ' << grmax[0] << endl;
	
	for (int i=0,j=0;i < int (esp.size()); ++i, ++j){
	  if (j==2) {opl << endl; j=0;}
	  opl << ' ' << esp[i] << ' ' ; 
	}
	
	
	
	
	opl.close(); 
	
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

