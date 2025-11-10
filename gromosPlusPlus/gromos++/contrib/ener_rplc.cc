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
 * @file ener_rplc.cc
 * recalculate interaction energies after replacing a soft-core
 * interaction site with a given molecule
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor ener_rplc
 * @section ener_rplc nonbonded energies after replacement of atoms
 * @author @ref co
 * @date 26.02.2007 
 *

 * calculate the time series of interaction energies after replacing a single
 * solute atom with a given molecule. This is a very specific program to 
 * calculate solvation free energies from single-step perturbation 
 * calculations as described in [J. Comput. Chem. 20 (1999) 1604-1617] and 
 * [J. Phys. Chem. B 105 (2001) 11264-11274]. It assumes a trajectory file from
 * a simulation of a single soft reference site in solvent. The user specifies
 * a topology with a solute molecule as well as coordinates for these. 
 * The centre of geometry of specified fit-atoms is placed on the soft 
 * reference site and the interaction energies of the solute atoms is 
 * calculated. Additional sampling of the solute-atoms can be performed in the
 * form of rotations and translations.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;topology&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms" for energy calculation&gt; </td></tr>
 * <tr><td> \@fitatoms</td><td>&lt;@ref AtomSpecifier "atoms" to fit to (cog)&gt; </td></tr>
 * <tr><td> \@trans</td><td>&lt;nr translations&gt; &lt;max distance&gt; </td></tr>
 * <tr><td> \@rot</td><td>&lt;nr rotations&gt; </td></tr>
 * <tr><td> \@seed</td><td>&lt;random number seed&gt; </td></tr>
 * <tr><td> \@solute</td><td>&lt;coordinates of the molecule to put in&gt; </td></tr>
 * <tr><td> \@cut</td><td>&lt;cut-off distance&gt; </td></tr>
 * <tr><td> \@eps</td><td>&lt;epsilon for reaction field correction&gt; </td></tr>
 * <tr><td> \@kap</td><td>&lt;kappa for reaction field correction&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  ener_rplc
    @topo
    @pbc
    @atoms
    @fitatoms
    @trans
    @rot
    @seed
    @solute
    @cut
    @eps
    @kap
    @traj
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <math.h>
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/Box.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/Energy.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Matrix.h"
#include "../src/gromos/Exception.h"

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;


void rotate(System &sys, AtomSpecifier as, vector<Vec> coord, Vec angle);
void translate(System &sys, AtomSpecifier as, vector<Vec> coord, Vec trans);

int main(int argc, char **argv){

  Argument_List knowns;
  knowns << "topo" << "pbc" << "atoms" << "fitatoms" << "trans" << "rot" 
         << "seed" << "solute" << "traj" << "cut" << "eps" << "kap";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo     <topology>\n";
  usage += "\t@pbc      <boundary type>\n";
  usage += "\t@atoms    <atoms for energy calculation>\n";
  usage += "\t@fitatoms <atoms to fit to (cog)>\n";
  usage += "\t@trans    <nr translations> <max distance>\n";
  usage += "\t@rot      <nr rotations>\n";
  usage += "\t@seed     <random number seed>\n";
  usage += "\t@solute   <coordinates of the molecule to put in>\n";
  usage += "\t@cut      <cut-off distance>\n";
  usage += "\t@eps      <epsilon for reaction field correction>\n";
  usage += "\t@kap      <kappa for reaction field correction>\n";
  usage += "\t@traj     <trajectory files>\n";
  
 
  try{
    Arguments args(argc, argv, knowns, usage);

    //  read topology, this is a topology with one solute molecule
    InTopology it(args["topo"]);
    System sysw(it.system());
    GromosForceField gff(it.forceField());

    // read it in
    InG96 ic(args["solute"]);
    ic >> sysw;
    ic.close();
    
    System refSys(sysw);

    // gather it
    Boundary *pbcw = BoundaryParser::boundary(sysw, args);
    // GatherParser
    Boundary::MemPtr gathmethodw = args::GatherParser::parse(sysw,refSys,args);
    (*pbcw.*gathmethodw)();

    // create a topology with the soft atom
    MoleculeTopology mt;
    AtomTopology at;
    mt.addAtom(at);
    System sys;
    sys.addMolecule(mt);
    sys.addSolvent(sysw.sol(0));
    
    // and a third system
    System sysnew(sysw);

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sysnew, args);
    // GatherParser
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

    // read in atoms to consider for the energy calculation
    utils::AtomSpecifier as(sysnew);
    {
      Arguments::const_iterator iter=args.lower_bound("atoms"), 
	to=args.upper_bound("atoms");
      for(;iter!=to; ++iter)
	as.addSpecifier(iter->second);
    }
 
    // read in atoms to consider for fitting
    utils::AtomSpecifier fit(sysw);
    {
      Arguments::const_iterator iter=args.lower_bound("fitatoms"),
	to=args.upper_bound("fitatoms");
      for(;iter!=to; ++iter)
	fit.addSpecifier(iter->second);
    }
    
    // define energy class
    Energy en(sysnew, gff, *pbc);
    
    // read in the energy things
    //   get cut-off distance

    en.setCutOff(args.getValue<double>("cut", false, 1.4));
    en.setRF(args.getValue<double>("eps", false, 1.0),
            args.getValue<double>("kap", false, 0.0));

    // tel energy which atoms to consider
    en.setAtoms(as);
    
    // read in number of translations and the maximum translation
    vector<double> transarg = args.getValues<double>("trans", 2, false,
            Arguments::Default<double>() << 1.0 << 0.05);
    int numtrans = int(transarg[0]);
    double rmax = transarg[1];

    
    // read in number of rotations
    int numrot = args.getValue<int>("rot", false, 1);

    // set the random number generator seed
    srand(args.getValue<int>("seed"));
    
    // count the number of frames
    int numframes=0;

    // write a title
    cout << "#" 
	 << setw(5) << "frm"
	 << setw(6) << "rot"
	 << setw(6) << "trans"
	 << setw(22) << "vdw"
	 << setw(22) << "crf"
	 << setw(22) << "tot"
	 << endl;
    
    // loop over all trajectories
    for(Arguments::const_iterator iter=args.lower_bound("traj"), 
	  to=args.upper_bound("traj"); iter!=to; ++iter){

      ic.open(iter->second);
      ic.select("ALL");
      
      // loop over single trajectory
      while(!ic.eof()){
	ic >> sys;

	numframes++;
	
	// sys now contains the coordinates of the soft atoms
	// to make things really slow, copy it over

	// move the solute to the soft cavity
	Vec cog(0.0,0.0,0.0);
	for(int i=0; i<fit.size(); ++i)
	  cog+=fit.pos(i);
	cog/=fit.size();
	
	for(int i=0; i<sysw.mol(0).numAtoms(); i++){
	  sysnew.mol(0).pos(i) = sysw.mol(0).pos(i) 
	    - cog
	    + sys.mol(0).pos(0);
	}
	//cout << sysnew.sol(0).numPos() << endl;
	//sysnew.sol(0).setnumCoords(sys.sol(0).numCoords());
	sysnew.sol(0).setNumPos(0);
	
	//cout << sysnew.sol(0).numCoords() << endl;
	for(int i=0; i<sys.sol(0).numPos(); i++){
	  //cout << i<< endl;
	  sysnew.sol(0).addPos(sys.sol(0).pos(i));
	  
	  //sysnew.sol(0).pos(i) = sys.sol(0).pos(i);
	}
	
	//cout << sysnew.sol(0).numPos() << endl;
	
	//System sysnew(sysw);
	//sysnew.addSolvent(sys.sol(0));
	sysnew.box().K() = sys.box().K();
        sysnew.box().L() = sys.box().L();
        sysnew.box().M() = sys.box().M();
        //for(int i=0; i< 3; i++) sysnew.box()[i] = sys.box()[i];
	sysnew.hasBox = true;
	
	// gather (with any method) to connect the molecules
	(*pbc.*gathmethod)();
	
	// store the original coordinates in a vector of Vec
	vector<Vec> coord;
	for(int k=0; k<as.size(); k++)
	  coord.push_back(as.pos(k));
	
	// define angles to rotate. The first rotation is zero
	// this is expensive, but we will hardly notice it
	Vec angle(0.0,0.0,0.0);
	int irot=0;
	
	do{
	  //every rotation starts from the original coordinates
	  rotate(sysnew, as, coord, angle);
		  
	  // store these rotated coordinates, so that we can translate
	  // from the same coordinates all the time
	  vector<Vec> coord2;
	  for(int k=0; k<as.size(); k++)
	    coord2.push_back(as.pos(k));
	  
	  //cout << "size " << coord2.size();
	  
	  // again we start with a zero translation
	  Vec trans(0.0,0.0,0.0);
	  int itrans=0;
	  
	  do{
	    translate(sysnew, as, coord2, trans);
	    
	    // now calculate the energies
	    
	    en.calc();
	    cout.setf(ios::right, ios::adjustfield);
	    cout << setw(6) << numframes
		 << setw(6) << irot + 1
		 << setw(6) << itrans + 1
		 << setw(22) << en.vdw()
		 << setw(22) << en.el()
		 << setw(22) << en.tot()
		 << endl;
	    
	    // and prepare then next translation
	    trans=Vec(double(rand()), double(rand()), double(rand()));
	    trans*=2.0*rmax/RAND_MAX;
	    trans-=Vec(rmax, rmax, rmax);
            itrans++;
	    
	  }while(itrans<numtrans);
	  
	  // and prepare the next rotation 
	  angle=Vec(double(rand()), double(rand()), double(rand()));
	  angle*=2.0*M_PI/RAND_MAX;
	  //cout << angle[0] << " " << angle[1] << " " << angle[2] << endl;
	  
	  irot++;
	  
	}while(irot<numrot);
      }
    }
    
    ic.close();
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
	
void translate(System &sys, AtomSpecifier as, vector<Vec> coord, Vec trans)
{
  //this will be a great function
  for(int i=0; i<as.size(); i++)
    sys.mol(as.mol(i)).pos(as.atom(i))=coord[i]+trans;
}

void rotate(System &sys, AtomSpecifier as, vector<Vec> coord, Vec angle)
{
  // this one will be a bit more complex

  // to pretend we care about speed, first calculate the cosines and sines
  // of all angles
  Vec angle_cos, angle_sin;
  
  for(int i=0; i<3; i++){
    angle_cos[i]=cos(angle[i]);
    angle_sin[i]=sin(angle[i]);
  }
  // now we copy things around
  // first calculate the centre of geometry
  Vec com(0.0,0.0,0.0);
  double mass, totmass=0.0;
  
  for(int k=0; k<as.size(); k++){
    mass=sys.mol(as.mol(k)).topology().atom(as.atom(k)).mass();
    com+=coord[k]*mass;
    totmass+=mass;
  }
  com/=totmass;
 
  // and move them to that point
  for(int k=0; k<as.size(); k++)
    coord[k]-=com;
  
  // The product of three rotation matrices about three axes
  /*
   * (   1.0   0.0   0.0)   ( cosy   0.0   siny)   ( cosx   sinx   0.0)
   * (   0.0  cosz  sinz) X (  0.0   1.0    0.0) X (-sinx   cosx   0.0)
   * (   0.0 -sinz  cosz)   (-siny   0.0   cosy)   (  0.0    0.0   1.0)  
   */

  gmath::Matrix rot(3,3);
  rot(0,0)=angle_cos[0]*angle_cos[1];
  rot(1,0)=-angle_sin[0]*angle_cos[2]-angle_cos[0]*angle_sin[1]*angle_sin[2];
  rot(2,0)=angle_sin[0]*angle_sin[2]-angle_cos[0]*angle_sin[1]*angle_cos[2];
  rot(0,1)=angle_sin[0]*angle_cos[1];
  rot(1,1)=angle_cos[0]*angle_cos[2]-angle_sin[0]*angle_sin[1]*angle_sin[2];
  rot(2,1)=-angle_cos[0]*angle_sin[2]-angle_sin[0]*angle_sin[1]*angle_cos[2];
  rot(0,2)=angle_sin[1];
  rot(1,2)=angle_cos[1]*angle_sin[2];
  rot(2,2)=angle_cos[1]*angle_cos[2];

  //rotate the coordinates and move them back to their original place
  for(int k=0; k<as.size(); k++)
    sys.mol(as.mol(k)).pos(as.atom(k))=rot*coord[k] + com;
}

