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
 * @file unify_box.cc
 * Convert boxshapes
 */

/**
 * @page programs Program Documentation
 *
 * @anchor unify_box
 * @section unify_box Convert boxshapes
 * @author @ref mc
 * @date 05-03-2003
 *
 * Program unify box can convert different boxshapes. All periodic boxes can be
 * described as a triclinic box, which is defined by vectors K,L and M. The
 * program is mostly used to convert a truncated octahedral box into a
 * triclinic box or vice versa, according to [H. Bekker, J Comp Chem, 18 (15),
 * 1930, 1997]. The user can also specify a rotation matrix and K,L and M 
 * vectors directly.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> [\@from_pbc</td><td>&lt;original boundary condition&gt;] </td></tr>
 * <tr><td> \@to_pbc</td><td>&lt;target boundary condition&gt; </td></tr>
 * <tr><td> [\@rot</td><td>&lt;rotation matrix&gt;] </td></tr>
 * <tr><td> [\@KLM</td><td>&lt;K-, L-, and M-vectors&gt;] </td></tr>
 * <tr><td> \@pos</td><td>&lt;coordinate file&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  unify_box
    @topo       ex.top
    @from_pbc   t
    @to_pbc     c
    @pos        exref.coo
 @endverbatim
 *
 * <hr>
 */

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Matrix.h"
#include "../src/gromos/Exception.h"

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace gmath;
using namespace std;

// should go to the matrix class
std::ostream &operator<<(std::ostream &os, const Matrix &m)
{
  os << "[\n";
  for(int i=0; i<3;++i){
    os << "[";
    for(int j=0; j<3; ++j){
      os << std::setw(20) << m(i,j);
    }
    os << std::setw(5) <<  "]\n";
  }
  os << "]\n";

  return os;
  
}
void rot_trunc_tric(gcore::Box & box, gmath::Matrix & rot, 
		    gcore::Box &newbox, bool forward);

int main(int argc, char **argv){

  Argument_List knowns;
  knowns << "topo" << "from_pbc" << "to_pbc" << "rot" << "KLM" << "pos";

  std::string usage = argv[0];
  usage += "\n\t@topo       <molecular topology file>\n";
  usage += "\t[@from_pbc <original boundary condition>]\n";
  usage += "\t@to_pbc    <target   boundary condition>\n";
  usage += "\t[@rot      <rotation matrix>]\n";
  usage += "\t[@KLM      <K-, L-, and M-vectors>]\n";
  usage += "\t@pos       <coordinate file>\n";
 
  try{
    Arguments args(argc, argv, knowns, usage);

    //  read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    
    Boundary *to_pbc = BoundaryParser::boundary(sys, args, "to_pbc");

    InG96 ic(args["pos"]);
    ic.select("ALL");
    ic >> sys;
    
    // output
    OutG96S oc(std::cout);
    oc.select("ALL");
    std::string title=ic.title();

    gmath::Matrix rot(3,3);
    gcore::Box newbox;
    
    if(args.count("from_pbc") > 0){
      if(args.count("rot")>=0 || args.count("KLM")>=0)
	throw gromos::Exception("unify_box",
				"Either give @rot and @KLM or @from_pbc, not both");
      std::string from_pbc=args["from_pbc"];
      if(from_pbc=="t"){
	if(to_pbc->type()=='c') {
	  rot_trunc_tric(sys.box(), rot, newbox, true);
	  title+="converted from truncated octahedron to triclinic box";
	}
	else
	  throw gromos::Exception("unify_box",
				  "Feature not implemented");
      }
      else if(from_pbc=="c"){
	if(to_pbc->type()=='t') {
	  rot_trunc_tric(sys.box(), rot, newbox, false);
	  title+="converted from triclinic box to truncated octahedron";
	}
	else
	  throw gromos::Exception("unify_box",
				  "Feature not implemented");
      }
      else
	throw gromos::Exception("unify_box",
				"Feature not implemented");
    }
    else{
      //we have to read in the rot and the KLM
      if(args.count("KLM")<=0)
	throw gromos::Exception("unify_box",
				"If you give rot, then you should give KLM");
      if(args.count("rot")!=9)
	throw gromos::Exception("unify_box",
				"A matrix has nine (9) elements");
      Arguments::const_iterator iter=args.lower_bound("rot");
      for(int i=0; i< 3; ++i){
	for(int j=0; j<3; ++j, ++iter){
	  rot(i,j)=atof(iter->second.c_str());
	}
      }
      //      std::cout << "the matrix we read\n" << rot << std::endl;
      
      if(args.count("KLM")!=9)
	throw gromos::Exception("unify_box",
				"We expect nine elements for three vectors"
				", but in which order?");
      iter=args.lower_bound("KLM");
      for(int i=0; i< 3; ++i){
	newbox.K()[i]=atof(iter->second.c_str());
	++iter;
	newbox.L()[i]=atof(iter->second.c_str());
	++iter;
	newbox.M()[i]=atof(iter->second.c_str());
	++iter;
      }
      newbox.update_triclinic();
      
      std::cerr << "Lattice vectors:\n"
	    << "K" << gmath::v2s(newbox.K())
	    << "\nL" << gmath::v2s(newbox.L())
	    << "\nM" << gmath::v2s(newbox.M())
	    << "\n";
      title+="converted on your own whim";
    }
    
    sys.box()=newbox;
    
    gmath::Vec origo(0.0,0.0,0.0);
    
    for(int i=0;i<sys.numMolecules();i++){
      for(int j=0;j<sys.mol(i).topology().numAtoms();j++){
	// rotate the coordinate system
	sys.mol(i).pos(j) = rot * sys.mol(i).pos(j);
	// and take the nearest image with respect to the origo
	sys.mol(i).pos(j) = to_pbc->nearestImage(origo, sys.mol(i).pos(j),
						 sys.box());
      }
    }
    for(int j=0; j<sys.sol(0).numPos(); j++){	
      // rotate the coordinate system
      sys.sol(0).pos(j) = rot * sys.sol(0).pos(j);
      // and take the nearest image with respect to the origo
      sys.sol(0).pos(j) = to_pbc->nearestImage(origo, sys.sol(0).pos(j),
					       sys.box());
    }
    
    if(to_pbc->type()=='c') {
      sys.box().setNtb(gcore::Box::triclinic);
      sys.box().boxformat()=gcore::Box::triclinicbox;
    }
    

    // write the new coordinates
    oc.writeTitle(title);
    oc << sys;
    
    ic.close();
    oc.close();
  }
  catch (const gromos::Exception &e){
    std::cerr << e.what() << std::endl;
    return 1;
  }
  return 0;
}

void rot_trunc_tric(gcore::Box & box, gmath::Matrix & rot, 
		    gcore::Box &newbox, bool forward)
{
  const double third = 1.0 / 3.0;
  const double sq3i = 1.0/sqrt(3.0);
  const double sq2i = 1.0/sqrt(2.0);

  if(forward){
    const double d = 0.5*sqrt(3.0) * box.K()[0];
    //const double d = 0.5*sqrt(3.0) * box[0];
    
    newbox.K()[0] =  d;
    newbox.K()[1] =  0.0;
    newbox.K()[2] =  0.0;
    
    newbox.L()[0] =  third * d;
    newbox.L()[1] =  2 * third * sqrt(2.0) * d;
    newbox.L()[2] =  0.0;
    
    newbox.M()[0] = -third * d;
    newbox.M()[1] =  third * sqrt(2.0) * d;
    newbox.M()[2] =  third * sqrt(6.0) * d;
    
    newbox.update_triclinic();
    std::cerr << "K_L_M       : " << newbox.K_L_M() << std::endl;
    std::cerr << "cross_K_L_M 1: " << gmath::v2s(newbox.cross_K_L_M()[0]) << std::endl;
    std::cerr << "cross_K_L_M 2: " << gmath::v2s(newbox.cross_K_L_M()[1]) << std::endl;
    std::cerr << "cross_K_L_M 3: " << gmath::v2s(newbox.cross_K_L_M()[2]) << std::endl;
  }
  else{
    const double d = 2*box.K()[0] / sqrt(3.0);
    
    newbox.K()[0]=d;
    newbox.K()[1]=0;
    newbox.K()[2]=0;

    newbox.L()[0]=0;
    newbox.L()[1]=d;
    newbox.L()[2]=0;

    newbox.M()[0]=0;
    newbox.M()[1]=0;
    newbox.M()[2]=d;
  }

  // taking K = box*(.5, .5, .5)
  //        L = box*(-.5, .5, .5)
  //        M = box*(-.5, -.5, .5)

 
  // now rotate the truncated octahedron !!!
  // the analytical rotation matrix
  
  // we need to get Ky = Kz = Lz = 0.0
  // (the other conditions are fulfilled automatically by truncoct) MAYBE

  rot = Matrix(Vec(sq3i, -2*sq2i*sq3i, 0),
	       Vec(sq3i, sq3i*sq2i, -sq2i),
	       Vec(sq3i, sq2i*sq3i, sq2i));
  

  // to get a triclinic cell back to a truncated octahedron
  // we need the inverse of the rotation matrix
  // this happens to be the transpose...
  if(forward){
    std::cerr << "Converting from truncated octahedron to triclinic box\n\n";
  }
  else{
    rot = rot.transpose();
    std::cerr << "Converting from  triclinic box to truncated octahedron\n\n";
  }
  std::cerr << "Lattice vectors:\n"
	    << "K" << gmath::v2s(newbox.K())
	    << "\nL" << gmath::v2s(newbox.L())
	    << "\nM" << gmath::v2s(newbox.M())
	    << "\n";
  
  std::cerr << "Rotation matrix\n" << rot << std::endl;
  
  
  // write out the diagonals
  // the cut-off has to be less than 1/2 * smallest diagonal
  std::cerr << "Diagonals:" << std::endl;
  
  Vec diag = newbox.K() + newbox.L() + newbox.M();
  std::cerr << "\t+++\t" << diag.abs() << std::endl;
  
  diag = newbox.K() + newbox.L() - newbox.M();
  std::cerr << "\t++-\t" << diag.abs() << std::endl;
  
  diag = newbox.K() - newbox.L() + newbox.M();
  std::cerr << "\t+-+\t" << diag.abs() << std::endl;
  
  diag = -newbox.K() + newbox.L() + newbox.M();
  std::cerr << "\t-++\t" << diag.abs() << std::endl;
}
