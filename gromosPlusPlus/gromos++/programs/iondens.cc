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
 * @file iondens.cc
 * monitor the average density of ions in the simulation box
 */

/**
 * @page programs Program Documentation
 *
 * @anchor iondens
 * @section iondens monitor the average density of ions in the simulation box
 * @author @ref mk
 * @date 23.6.06
 *
 * Program iondens calculates the average density of ions (or other particles)
 * over a trajectory file. A rotational fit of the system onto the solute
 * can be performed to correct for rotations of the complete simulation box. The
 * density will be calculated on a grid of points. Atoms specified by \@atoms are 
 * use to center the grid and for the rotational fit. If \@atoms is not specified, 
 * no rotational fit will be performed, the centre of the grid will be calculated 
 * by the centre of geometry of 1:a of the system. Two sets of densities can be
 * written out, containing 1) occupancies on the grid points, relative to the
 * maximally occupied gridpoint (in file grid.pdb), or 2) occupancies as a 
 * percentage of the number of frames (in file gridnf.pdb). These two output files 
 * are written in pdb format. User specified cutoffs determine which gridpoints  
 * will be written out.
 * 
 * The program iondens also writes out the average structure of the system over
 * the trajectories in pdb format in a file "aver.pdb".
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@grspace</td><td>&lt;grid spacing (default: 0.2 nm)&gt; </td></tr>
 * <tr><td> \@ions</td><td>&lt;@ref AtomSpecifier "ions" to monitor&gt; </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms" to use for fit&gt; </td></tr>
 * <tr><td> \@ref</td><td>&lt;reference coordinates&gt; </td></tr>
 * <tr><td> \@thresholds</td><td>&lt;threshold values for occupancy percentages (default: 20 and 5)&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  iondens
    @topo       ex.top
    @pbc        r
    @grspace    0.2
    @ions       2-44:1
    @atoms      1:CA
    @ref        exref.coo
    @thresholds 20 5
    @traj       ex.tr
 @endverbatim
 *
 * <hr>
 */
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/args/ReferenceParser.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/gcore/Box.h"
#include "../src/gio/OutPdb.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/gromos/Exception.h"

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace fit;
using namespace std;
using namespace gmath;


int main(int argc, char **argv){

  Argument_List knowns;
  knowns << "topo" << "pbc" << "atoms" << "ref" << "grspace" << "traj" 
         << "ions" << "thresholds";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo         <molecular topology file>\n";
  usage += "\t@pbc          <boundary type>\n";
  usage += "\t[@grspace      <grid spacing (default: 0.2 nm)>]\n";
  usage += "\t@ions         <ions to monitor>\n";
  usage += "\t[@atoms       <atoms to use for fit>]\n";
  usage += "\t[@ref          <reference coordinates>]\n";
  usage += "\t@thresholds   <threshold values for occupancy percentages (default: 20 and 5)>\n";
  usage += "\t@traj         <trajectory files>\n";
  
 
  try{
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    // make two systems
    System sys(it.system());
    System refSys(it.system());

    // which atoms do we use to center the grid on and for the rotational fit?
    utils::AtomSpecifier atomsfit(refSys);
    utils::AtomSpecifier atoms(sys);
    bool dofit=false;
    for(Arguments::const_iterator it=args.lower_bound("atoms");
	it!= args.upper_bound("atoms"); ++it){
      atomsfit.addSpecifier(it->second);
    }
    if(atomsfit.size()==0){
      dofit=false;
      cout << "# No rotational fit performed" << endl;
      atoms.addSpecifier("1:a");
    }
    else{
      dofit=true;
      cout << "# Performing a rotational fit" << endl;
      
      // to get the same atomlist atomsift -> atoms 
      //atoms=atomsfit;
      const vector<string> & spec = atomsfit.toString();

      for(vector<string>::const_iterator it = spec.begin(), to = spec.end();
              it != to; ++it) {
        atoms.addSpecifier(*it);    
      }
      //atoms.setSystem(sys);
    }

    //which ions to consider
    utils::AtomSpecifier ions(sys);
    for(Arguments::const_iterator it=args.lower_bound("ions");
	it!=args.upper_bound("ions"); ++it){
      ions.addSpecifier(it->second);
    }
    
    if(ions.size()==0){
      throw gromos::Exception("iondens", "No ions specified to monitor");
    }
    
    // get grid spacing
    double space = args.getValue<double>("grspace", false, 0.2);
    // get threshold values
    vector<double> thres = args.getValues<double>("thresholds", 2, false,
            Arguments::Default<double>() << 20.0 << 5.0);

    // input coordinate
    InG96 ic;
    
    if(args.count("ref")==1)
      ic.open(args["ref"]);
    else
      ic.open(args["traj"]);
    ic.select("ALL");   
    ic >> refSys;
    ic.close();    
    
    cout << "# Monitoring the positions of " << ions.size() << " ions" << endl;

    // parse boundary conditions for the reference system
    Boundary *pbc = BoundaryParser::boundary(refSys, args);

    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

    (*pbc.*gathmethod)();

    delete pbc;

    // now set the pbc for the system
    pbc = BoundaryParser::boundary(sys, args);

    Reference reffit(&refSys);
    reffit.addAtomSpecifier(atomsfit);

    RotationalFit rf(&reffit);

    Reference ref(&sys);
    ref.addAtomSpecifier(atoms);
    
    // wirte out reference structure, we don't need this here.    
    //OutCoordinates *oc;
    //oc = new OutPdb();
    //ofstream os("ref.pdb");
    //OutPdb opdb(os);
    //opdb << refSys;
    //os.close();

    // create a system for the average structure, set the coords to 0.0
    System aver(sys);
    for (int i=0;i < aver.numMolecules(); ++i){
      for (int j=0; j < aver.mol(i).numAtoms(); ++j){
	aver.mol(i).pos(j)=Vec(0.0, 0.0, 0.0);
      }
    }

    // set some vectors and stuff that we may need.
    Vec center (0.0, 0.0, 0.0), cog (0.0,0.0,0.0), bx(0.0,0.0,0.0),
      grmin(0.0,0.0,0.0), grmax(0.0,0.0,0.0);

    int nx=0,ny=0,nz=0, numFrames=0;

    // build grid positions
    vector<Vec> densgrid;

    Box box = refSys.box();
    box.stretch_K(1.2 * box.K().abs());
    box.stretch_L(1.2 * box.L().abs());
    box.stretch_M(1.2 * box.M().abs());

    nx=int(rint(box.K().abs()/space)); ny=int(rint(box.L().abs()/space)); nz=int(rint(box.M().abs()/space));
    //nx=int(rint(box[0]/space)); ny=int(rint(box[1]/space)); nz=int(rint(box[2]/space));

    //make sure we have equal number of points in x,y,z
    //dont know if this is nescessary for the programs that read in grids
    //i guess it wont hurt either...
    nx=max(nx,ny); nx=max(nx,nz);
    ny=nx; nz=ny;

    bx[0]=nx*space; bx[1]=ny*space; bx[2]=nz*space;

    Vec start=-bx/2;
 
    for(int i=0;i<nx;i++){
      for(int j=0;j<ny;j++){
	for(int k=0;k<nz;k++){ 
	  Vec tmp;      
	  tmp[0]=start[0]+space*i;
	  tmp[1]=start[1]+space*j;
	  tmp[2]=start[2]+space*k;
	  densgrid.push_back(tmp);
	}
      }
    }
    grmin=start;
    grmax=densgrid[densgrid.size()-1];

    vector<int> ioncount; 
    
    for (int i=0; i < int (densgrid.size()); ++i){
        ioncount.push_back(0);
    }
     
    // loop over all trajectories
    for(Arguments::const_iterator 
	  iter=args.lower_bound("traj"),
	  to=args.upper_bound("traj");
	iter!=to; ++iter){

      // open file
      ic.open((iter->second).c_str());
      ic.select("ALL");
      // loop over single trajectory
      while(!ic.eof()){
	numFrames+=1;
	
	ic >> sys;	

	(*pbc.*gathmethod)();
	if(dofit)
	  rf.fit(&sys);	   
	
        // sum average position
	for (int i=0;i < aver.numMolecules(); ++i){
	  for (int j=0; j < aver.mol(i).numAtoms(); ++j){
	    aver.mol(i).pos(j)+=sys.mol(i).pos(j);
	  }
	}
        
        // calculate cog again for the selected atoms
	cog=PositionUtils::cog(sys, ref);
        
	// calculate nim with respect to cog from above
	for (unsigned int i=0; i < ions.size(); ++i){
	  ions.pos(i) = pbc->nearestImage(cog,ions.pos(i),sys.box());
	}
	
	double rmin = 1000000, r=0; int p=0;
	// determine ion position with respect to grid
	for (int i=0; i< int(ions.size()); ++i){
	  p=0, rmin = 1000000;
	  Vec ion = ions.pos(i);
	  if(fabs(ion[0]-start[0])>bx[0] || 
	     fabs(ion[1]-start[1])>bx[1] ||
	     fabs(ion[2]-start[2])>bx[2])
	    cout << "#Ion found outside the grid!\n\t" 
		 << ion[0] << " " << ion[1] << " " << ion[2] << endl;
	  
	  for (int j=0; j < int (densgrid.size()); ++j){
	    Vec tmp = densgrid[j];
	    r = (tmp-ion).abs();
	    if (rmin > r) {
	      rmin = r;
	      p = j;}
	  }
	  ioncount[p]+=1;
	}

      }
      ic.close();
    }
    
    // average
    for (int i=0;i < aver.numMolecules(); ++i){
      for (int j=0; j < aver.mol(i).numAtoms(); ++j){
	aver.mol(i).pos(j)=aver.mol(i).pos(j)/numFrames;
      }
    }
    // write out average structure
    ofstream oss("aver.pdb");
    OutPdb oopdb(oss);
    oopdb << aver;
    oss.close();

    // get max_element in ioncount
    vector<int>::iterator maxel = std::max_element(ioncount.begin(),ioncount.end());
  
    // normalize * 100
    vector<double> per;
    vector<double> pernf;
    for (int i=0; i < int (ioncount.size()); ++i){ 
      per.push_back(100.0 * double(ioncount[i])/double(*maxel));
      pernf.push_back(100.0 * double(ioncount[i])/double(numFrames*ions.size()));
      
      //do something stupid...  WHY??
      //if (per[i] == 100) {per[i]=99.99;}
      //if (pernf[i] == 100) {pernf[i]=99.99;}
    }

    //write out 'da grid with b-facs
    ofstream perfile;

    perfile.open("grid.pdb");
    int count=0;
    for (int i=0; i < int (ioncount.size()); ++i){ 
      if (per[i] > thres[0]) {
	++count;
	Vec tmp = densgrid[i];
	perfile.setf(ios::fixed, ios::floatfield);
	perfile.setf(ios::unitbuf);
	perfile.precision(3);

	perfile << "ATOM";
	perfile.setf(ios::right, ios::adjustfield);
	perfile << setw(7) << count;
	perfile.setf(ios::left, ios::adjustfield);
	perfile << "  " <<setw(4) << "AR";
	perfile << setw(4) << "AR";
	perfile.setf(ios::right, ios::adjustfield);
	perfile << setw(5) << "1" << "    "
		<< setw(8) << tmp[0]*10
		<< setw(8) << tmp[1]*10
		<< setw(8) << tmp[2]*10
		<< "  1.00 "; 
	perfile.precision(2);
	perfile << per[i]
		<< endl;
      }
    }
/*
    perfile.open("grid.dat");
    perfile << "# (x, y, z) with respect to the cog of specified atoms"
            << "     "
            << "occupancy, relative to the maximally occupied gridpoint [%]"
            << endl;
    
    for (int i=0; i < int (ioncount.size()); ++i) {
      if (per[i] > thres[0]) {
        Vec tmp = densgrid[i];
        perfile.setf(ios::fixed, ios::floatfield);
        perfile.setf(ios::unitbuf);
        perfile.precision(3);
        perfile << setw(15) << tmp[0]
	        << setw(15) << tmp[1]
		<< setw(15) << tmp[2];
            
        perfile.precision(2);
        perfile << setw(20) << per[i]
	        << endl;
      }
    }
*/    
    perfile.close();

    //write out 'da grid with b-facs
    ofstream pernffile; 

    pernffile.open("gridnf.pdb");
    count=0;
    for (int i=0; i < int (ioncount.size()); ++i){ 
      if (pernf[i] > thres[1]) {
	++count;
	Vec tmp = densgrid[i];
	pernffile.setf(ios::fixed, ios::floatfield);
	pernffile.setf(ios::unitbuf);
	pernffile.precision(3);

	pernffile << "ATOM";
	pernffile.setf(ios::right, ios::adjustfield);
	pernffile << setw(7) << count;
	pernffile.setf(ios::left, ios::adjustfield);
	pernffile << "  " <<setw(4) << "AR";
	pernffile << setw(4) << "AR";
	pernffile.setf(ios::right, ios::adjustfield);
	pernffile << setw(5) << "1" << "    "
		  << setw(8) << tmp[0]*10
		  << setw(8) << tmp[1]*10
		  << setw(8) << tmp[2]*10
		  << "  1.00 "; 
	pernffile.precision(2);
	pernffile << pernf[i]
		  << endl;
      }
    }
/*
    pernffile.open("gridnf.dat");
    pernffile << "# (x, y, z) with respect to the cog of specified atoms"
              << "     "
              << "occupancy as a percentage of the number of frames [%]"
              << endl;
    
    for (int i=0; i < int (ioncount.size()); ++i){ 
      if (pernf[i] > thres[1]) {
        Vec tmp = densgrid[i];
        pernffile.setf(ios::fixed, ios::floatfield);
        pernffile.setf(ios::unitbuf);
        pernffile.precision(3);
        pernffile << setw(15) << tmp[0]
                  << setw(15) << tmp[1]
		  << setw(15) << tmp[2];
            
	pernffile.precision(2);
	pernffile << setw(20) << pernf[i]
	  	  << endl;
      }
    }
*/    
    pernffile.close();

  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

