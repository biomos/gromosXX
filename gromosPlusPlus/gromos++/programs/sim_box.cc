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
 * @file sim_box.cc
 * solvates a solute in a box of pre-equilibrated solvent
 */

/**
 * @page programs Program Documentation
 *
 * @anchor sim_box
 * @section sim_box solvates a solute in a box of pre-equilibrated solvent
 * @author @ref co
 * @date 7-6-07
 *
 * When simulating a molecule in solution or in a crystal containing solvent 
 * molecules, the atomic coordinates of the solvent molecules are to be
 * generated, if they are not available from experiment. Program sim_box can 
 * solvate a solute in a pre-equilibrated box of solvent molecules. The file
 * specifying the solvent configuration should contain a BOX block with the
 * dimensions corresponding to the pre-equilibrated density. The solvent 
 * topology is read from the solvent block in the specified topology.
 *
 * To prevent overlap between solute and solvent molecules, only solvent
 * molecules for which the centre of geometry is at a minimum distance from any
 * solute atom (which can be defined via the \@thresh flag) are put into the 
 * box. Before solvating the solute molecules, the solute can be rotated such
 * that the largest distance between any two solute atoms is directed along the
 * z-axes, and the largest atom-atom distance in the xy-plane lies in the
 * y-direction, by giving the \@rotate flag. The resulting box containing solute
 * and solvent molecules can be either rectangular or a truncated octahedron. 
 * Its dimensions can be specified directly via the \@boxsize flag. If this flag
 * is given, the box dimensions are read in from the BOX block in the solute
 * coordinate file. Alternatively, when the \@minwall flag is given, the solute
 * is put into a box filled with solvent molecules with box dimensions
 * guaranteeing a minimum distance between any solute molecule and the box
 * edges. Either one value can be specified for the \@minwall flag, resulting in
 * a cubic or truncated octahedron box, or three values can be specified, to
 * generate a rectangular box. In the latter case, the solute molecules can be 
 * gathered on request (by specifying the \@gather flag) and the \@rotate flag must be given.
 * 
 * Note that to solvate a solute in a triclinic box, one can use sim_box to 
 * generate a rectangular box and subsequently apply the appropriate symmetry
 * transformations on the generated box using the program @ref cry (see section 
 * V-2.17) or use sim_box to generate a truncated octahedral box and convert
 * that to a triclinic box using the program @ref unify_box (see section V-5.6).
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file of the solute&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;periodic boundary conditions (r or t)&gt; </td></tr>
 * <tr><td> \@pos</td><td>&lt;input coordinate file for the solute&gt; </td></tr>
 * <tr><td> \@solvent</td><td>&lt;input coordinate file for the pre-equilibrated solvent&gt; </td></tr>
 * <tr><td> [\@minwall</td><td>&lt;minimum solute to wall distance&gt;] </td></tr>
 * <tr><td> [\@thresh</td><td>&lt;minimum solvent to solute distance (default 0.23 nm)&gt;] </td></tr>
 * <tr><td> [\@boxsize</td><td>(use boxsize specified in input file)] </td></tr>
 * <tr><td> [\@gather</td><td>(gather solute)] </td></tr>
 * <tr><td> [\@rotate</td><td>(rotate solute: biggest axis along z, second along y)] </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  sim_box
    @topo     ex.top
    @pbc      r
    @pos      exref.co
    @solvent  h2o.g96
    @minwall  1.4
    @thresh   0.23
 @endverbatim
 *
 * <hr>
 */
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/Box.h"
#include "../src/bound/Boundary.h"
#include "../src/gio/InTopology.h"
#include "../src/fit/PositionUtils.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Matrix.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;
using namespace gmath;
using namespace bound;
using namespace utils;


void rotate_solute(System &sys, vector<double> &max, AtomSpecifier &as);
double calc_max_size(System &sys, AtomSpecifier &as, int dim);

int main(int argc, char **argv){

  Argument_List knowns;
  knowns << "topo" << "pbc" << "pos" << "solvent" << "minwall" << "thresh"
         << "boxsize" << "gather" << "rotate";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo      <molecular topology file>\n";
  usage += "\t@pbc       <periodic boundary conditions (r or t)>\n";
  usage += "\t@pos       <input coordinate file for the solute>\n";
  usage += "\t@solvent   <input coordinate file for the solvent>\n";
  usage += "\t[@minwall  <minimum solute to wall distance>]\n";
  usage += "\t[@thresh   <minimum solvent-solute distance (default 0.23 nm)>]\n";
  usage += "\t[@boxsize  (use boxsize specified in solute coordinate file)]\n";
  usage += "\t[@gather   (gather solute)]\n";
  usage += "\t[@rotate   (rotate solute: biggest axis along z, second along y)]\n";


  try{
    Arguments args(argc, argv, knowns, usage);
    
    // read topology make 2 systems 
    InTopology it(args["topo"]);
    System solu(it.system());

    System solv;
    solv.addSolvent(solu.sol(0));
    
    // set the definition of a hydrogen, for the heavy atom criterion
    for(int m=0; m<solu.numMolecules(); m++){
      solu.mol(m).topology().setHmass(1.008);
    }
    solu.sol(0).topology().setHmass(1.008);
    
    // read the minimum solute to wall distances.
    // three possibilities:
    // 1. nothing specified: the user will specify the box size
    //    using the @boxsize flag
    // 2. one number specified: the box will be cubic (pbc: r) or a
    //    trunc. oct. (pbc: t) with the size defined by maximum solute
    //    atom-atom distance plus twice this number
    // 3. three numbers specified: the box will be rectangular.
    //    the specified numbers are added to the
    //    maximum distances in x, y and z (after possibly rotating the solute)

    vector<double> minwall;
    {
      Arguments::const_iterator iter=args.lower_bound("minwall"), 
	    to=args.upper_bound("minwall");
      while(iter!=to){
	    minwall.push_back(atof(iter->second.c_str()));
	    ++iter;
      }
    }

    // read the minimum solvent-solute distance
    double minsol = args.getValue<double>("thresh", false, 0.23);
    double minsol2 = minsol * minsol;

    // check for the boxsize flag
    // if it is given, the box from the solute coordinates is used
    bool boxsize = false;
    if (args.count("boxsize") != -1)
      boxsize = true;
    
    // read coordinates into the systems
    InG96 ic;
    ic.open(args["pos"]);
    // we also read in any solvent that is already in the file
    ic.select("ALL");
    ic >> solu;
    ic.close();

    System refSys(it.system());

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(solu, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(solu,refSys,args);
    
    if (args.count("gather") != -1){
      // should not be used if boxsize is specified
      // ('cause you probably changed the boxsize by hand)
      if (boxsize)
	throw gromos::Exception("sim_box",
				"don't gather if you specified a boxsize in the solute coordinates\n"
				"they should be gathered before!");
      if (!solu.hasBox)
	throw gromos::Exception("sim_box",
				"solute gathering requested, but solute does not contain a box block!");
      (*pbc.*gathmethod)();
    }

    vector<double> max_dist;
    AtomSpecifier as(solu);

    if (args.count("rotate") != -1){
      rotate_solute(solu, max_dist, as);
    }
    else{
      if (minwall.size() == 3)
	throw gromos::Exception("sim_box",
				"solute should be rotated in order to align largest extension with z axis");
      
      max_dist.push_back(calc_max_size(solu, as, 3));
    }

    // determine box shape, calculate relevant things and check the 
    // input for consistency
    enum boundary_enum { vacuum, rectangular, cubic, truncoct, triclinic } boundary = vacuum;

    double size_corr=1.0;

    switch(pbc->type()){
      case('t'):
	boundary = truncoct;
	size_corr = 2.0 * sqrt(3.0)/3.0;
	
	if(minwall.size()>1)
	  throw(gromos::Exception("sim_box", 
	     "For truncated octahedral boxes you can only specify one number for @minwall"));

	if(minwall.size()==0){

          if (solu.box().K().abs() != solu.box().L().abs() || solu.box().K().abs() != solu.box().M().abs())
	  //if(solu.box()[0] != solu.box()[1] || solu.box()[0] != solu.box()[2])
	    throw(gromos::Exception("sim_box", 
				    "For truncated octahedral boxes, the specified boxsize should be the same in all dimensions"));
	}
	break;

      case('r'):
	boundary = rectangular;
	
	if (minwall.size() == 1)
          boundary = cubic;
        else if (minwall.size() == 0 && (solu.box().ntb() == gcore::Box::rectangular)
                && (solu.box().K().abs() == solu.box().L().abs()) && (solu.box().K().abs() == solu.box().M().abs()))
          boundary = cubic;

	break;

      case('c'):
	if(!boxsize)
	  throw gromos::Exception("sim_box",
       			  "boxsize has to be specified for triclinic");
	boundary = triclinic;
	
	break;
	
      case('v'):
	throw(gromos::Exception("sim_box", 
           "Why are you running this program if @pbc is vacuum?"));
	break;
    }
    
    // size of the solvent box
    vector<double> solvent_box(3, 0.0);
    
    if (boxsize){
      if (minwall.size())
	throw(gromos::Exception("sim_box", 
				"Cannot specify both boxsize and minwall."));

      if(!solu.hasBox)
	throw gromos::Exception("sim_box",
				"If you specify boxsize, the box dimensions should be in "
				"the BOX block of the solute");

      if(pbc->type() == 't'){
	if(solu.box().ntb() != gcore::Box::truncoct)
	  throw gromos::Exception("sim_box",
				  "Specified truncated octahedral for @pbc, and try to read boxsize from file\nbut no truncated octahedral specified there");
	
        solvent_box[0] = solu.box().K()[0];
        solvent_box[1] = solu.box().L()[1];
        solvent_box[2] = solu.box().M()[2];
        //for(int i=0; i<3; ++i){
	//  solvent_box[i] = solu.box()[0];
	//}
      }
      else if (pbc->type() == 'r'){
	if(solu.box().ntb() != gcore::Box::rectangular)
	  throw gromos::Exception("sim_box",
				  "Specified rectangular box for @pbc, and try to read boxsize from file\nbut no rectangular box specified there");
	

        solvent_box[0] = solu.box().K()[0];
        solvent_box[1] = solu.box().L()[1];
        solvent_box[2] = solu.box().M()[2];
        //for(int i=0; i<3; ++i){
	//  solvent_box[i] = solu.box()[i];
	//}
      }
      else if(pbc->type() == 'c'){
	if(solu.box().ntb() != gcore::Box::triclinic)
	  throw gromos::Exception("sim_box",
				  "Specified triclinic box for @pbc, and try to read boxsize from file\nbut no triclinic box specified there");

	// calculate the dimension of a large enough box to encompass the triclinic box
	// first construct the four diagonals
	vector<Vec> diag(4);
	diag[0] = 0.5*(solu.box().K() + solu.box().L() + solu.box().M());
	diag[1] = 0.5*(solu.box().K() + solu.box().L() - solu.box().M());
	diag[2] = 0.5*(solu.box().K() - solu.box().L() + solu.box().M());
	diag[3] = 0.5*(solu.box().K() - solu.box().L() - solu.box().M());
	
	// find the maximum and minimum x,y,z, store these in boxsize
	for(int j=0; j<3; ++j){
	  double maxdim=diag[0][j], mindim=diag[0][j];
	  for(int i=1; i<4; ++i){
	    if(diag[i][j] > maxdim) maxdim=diag[i][j];
	    if(diag[i][j] < mindim) mindim=diag[i][j];
	    if(-diag[i][j] > maxdim) maxdim= -diag[i][j];
	    if(-diag[i][j] < mindim) mindim= -diag[i][j];
	  }
	  solvent_box[j] = maxdim - mindim;
	}
      } // triclinic
      else{
	throw gromos::Exception("sim_box",
				"unknown boundary condition");
      }
      // check if all atoms are actually within the box.
      // If not, write a warning and put them in.
      bool warn_outside_box=false;
      
      for(int m=0; m < solu.numMolecules(); m++){
	for(int a=0; a < solu.mol(m).numAtoms(); a++){
	  
	  // are we inside the box
	  Vec check = pbc->nearestImage(Vec(0.0,0.0,0.0),
					solu.mol(m).pos(a) , solu.box());
	  
	  if(check[0]!=solu.mol(m).pos(a)[0] ||
	     check[1]!=solu.mol(m).pos(a)[1] || 
	     check[2]!=solu.mol(m).pos(a)[2]){
	    solu.mol(m).pos(a)=check;
	    warn_outside_box=true;
	    
	  }
	}
      }
      // and any waters that we already have based on the first atom
      for(int s=0; s < solu.numSolvents(); s++){
	for(int a=0; a < solu.sol(s).numPos(); 
	    a+=solu.sol(s).topology().numAtoms()){
	  // are we inside the box
	  Vec check = pbc->nearestImage(Vec(0.0,0.0,0.0),
					solu.sol(s).pos(a), solu.box());
	  
	  if(check[0]!=solu.sol(s).pos(a)[0] ||
	     check[1]!=solu.sol(s).pos(a)[1] || 
	     check[2]!=solu.sol(s).pos(a)[2]){
	    Vec shift=check-solu.sol(s).pos(a);
	    
	    for(int i=0; i< solu.sol(s).topology().numAtoms(); i++){
	      solu.sol(s).pos(a+i)+= shift;
	      warn_outside_box=true;
	    }
	  }
	}
      }

      if(warn_outside_box)
	cerr << "WARNING: not all atoms were within the specified box !\n"
	     << "         placed the atoms in the box before solvation\n"
	     << "         according to the specified box.\n";
      
    } // boxsize
    else{
      if (minwall.size() == 0)
	throw gromos::Exception("sim_box",
				"either use a specified boxsize "
				"or give a minimum distance from solute to the walls");
      
      if (boundary == truncoct){
	solu.box().setNtb(gcore::Box::truncoct);
	solu.box().K()[0] = size_corr*(max_dist[0] + 2 * minwall[0]);
        solvent_box[0] = solu.box().K()[0];
        solu.box().L()[1] = size_corr*(max_dist[0] + 2 * minwall[0]);
        solvent_box[1] = solu.box().L()[1];
        solu.box().M()[2] = size_corr*(max_dist[0] + 2 * minwall[0]);
        solvent_box[2] = solu.box().M()[2];
        //for(int i=0; i<3; i++){
	//  solu.box()[i] = size_corr*(max_dist[0] + 2 * minwall[0]);
	//  solvent_box[i] = solu.box()[i];
	//}
      }
      else{
	solu.box().setNtb(gcore::Box::rectangular);
	if (minwall.size() == 1){
          solu.box().K()[0] = max_dist[0]+2*minwall[0];
          solvent_box[0] = solu.box().K()[0];
          solu.box().L()[1] = max_dist[0]+2*minwall[0];
          solvent_box[1] = solu.box().L()[1];
          solu.box().M()[2] = max_dist[0]+2*minwall[0];
          solvent_box[2] = solu.box().M()[2];
	  //for(int i=0; i<3; i++){
	  //  solu.box()[i] = max_dist[0]+2*minwall[0];
	  //  solvent_box[i] = solu.box()[i];
	  //}
	}
	else{
	  // the solute has been rotated, 3 max_dist are known, 3 minwall distances are given
          solu.box().K()[0] = max_dist[2]+2*minwall[0];
	  solvent_box[0] = solu.box().K()[0];
          solu.box().L()[1] = max_dist[1]+2*minwall[1];
	  solvent_box[1] = solu.box().L()[1];
          solu.box().M()[2] = max_dist[0]+2*minwall[2];
	  solvent_box[2] = solu.box().M()[2];
	  //for(int i=0; i<3; i++){
	  //  solu.box()[i] = max_dist[2-i]+2*minwall[i];
	  //  solvent_box[i] = solu.box()[i];
	  //}
	}
      }
    }
    
    // read in the solvent coordinates. 
    // to make absolutely sure that there is a box block, check this    
    ic.open(args["solvent"]);
    ic.select("SOLVENT");
    ic >> solv;
    ic.close();

    if(!solv.hasBox)
      throw gromos::Exception("sim_box", 
			      "Could not read BOX block from solvent "
			      "coordinates");
    
    int num_solv_atoms_per_box = solv.sol(0).numPos();
	 
    // get the number of original solvent atoms    
    int numSolventAtoms=solu.sol(0).numPos();
      
    // move the solute to the centre of geometry
    Vec shiftcog=fit::PositionUtils::shiftToCog(&solu);
    // shift the solvent atoms in there as well
    // Is done in shiftcog nowadays!!
    //for(int i=0; i< numSolventAtoms; i++){
    //  solu.sol(0).pos(i)+= shiftcog;
    //}
    
    // calculate the solvent cog
    Vec solv_cog(0.0,0.0,0.0);
    for(int i=0; i<solv.sol(0).numPos(); i++)
      solv_cog+=solv.sol(0).pos(i);
    solv_cog/=solv.sol(0).numPos();
    // move to solvent cog
    for(int i=0; i<solv.sol(0).numPos(); i++)
      solv.sol(0).pos(i) -= solv_cog;

    // set the solute box size, calculate how many solvent boxes we should
    // use; calculate where to move the initial box
    vector<int> needed_boxes;
    Vec move_solvent;
    needed_boxes.push_back(int(solvent_box[0] / solv.box().K()[0]) + 1);
    move_solvent[0]= -0.5 * (needed_boxes[0] - 1) * solv.box().K()[0];
    needed_boxes.push_back(int(solvent_box[1] / solv.box().L()[1]) + 1);
    move_solvent[1]= -0.5 * (needed_boxes[1] - 1) * solv.box().L()[1];
    needed_boxes.push_back(int(solvent_box[2] / solv.box().M()[2]) + 1);
    move_solvent[2]= -0.5 * (needed_boxes[2] - 1) * solv.box().M()[2];
    //for(int i=0; i<3; i++){
    //  needed_boxes.push_back(int(solvent_box[i] / solv.box()[i]) + 1);
    //  move_solvent[i]= -0.5 * (needed_boxes[i] - 1) * solv.box()[i];
    //}

    // move the initial box so that after multiplication the complete box
    // is centered around the origin.
    for(int i=0; i< num_solv_atoms_per_box; i++)
      solv.sol(0).pos(i) += move_solvent;
    
    // do the multiplications
    for(int ix=0; ix < needed_boxes[0]; ix++){
      
      for(int iy=0; iy < needed_boxes[1]; iy++){
	
	for(int iz=0; iz < needed_boxes[2]; iz++){

	  if(ix != 0 || iy != 0 || iz != 0){
            Vec shift(ix * solv.box().K()[0], iy * solv.box().L()[1], iz * solv.box().M()[2]);
	    //Vec shift(ix * solv.box()[0], iy * solv.box()[1], iz * solv.box()[2]);

	    for(int atom = 0; atom < num_solv_atoms_per_box; atom++){
	      solv.sol(0).addPos(solv.sol(0).pos(atom) + shift);
	    }
	  }
	}
      }
    }

    int num_atoms_per_solvent=solv.sol(0).topology().numAtoms();
    int num_solvent_molecules=solv.sol(0).numPos() / num_atoms_per_solvent;

    // now we have to keep only those waters that are inside the box and 
    // far enough away from the solute
    // we look at the centre of geometry of the solvent molecule

    Vec o(0.0,0.0,0.0);
    double min_init =
      solvent_box[0] * solvent_box[0] +
      solvent_box[1] * solvent_box[1] +
      solvent_box[2] * solvent_box[2];
    
    for(int i=0; i< num_solvent_molecules; i++){

      // calculate the centre of geometry of this solvent
      Vec sol_i(0.0,0.0,0.0);
      for(int j=0; j< num_atoms_per_solvent; j++)
	sol_i += solv.sol(0).pos(num_atoms_per_solvent * i + j);
      sol_i /= num_atoms_per_solvent;
      
      // are we inside the box
      Vec check = pbc->nearestImage(o, sol_i, solu.box());

      if(check[0]==sol_i[0] && 
	 check[1]==sol_i[1] && 
	 check[2]==sol_i[2]){
	// yes we are in the box
	// calculate the closest distance to any solute
	double min2 = min_init;
	for(int m=0; m < solu.numMolecules(); m++){
	  for(int a=0; a < solu.mol(m).numAtoms(); a++){
	    if(! solu.mol(m).topology().atom(a).isH() && 
	       (check - solu.mol(m).pos(a)).abs2() < min2)
	      min2 = (check - solu.mol(m).pos(a)).abs2();
	  }
	}
	// or original solvent
	for(int j=0; j<numSolventAtoms; j++){
	  if(!solu.sol(0).topology().atom(j%num_atoms_per_solvent).isH() &&
	     (check - solu.sol(0).pos(j)).abs() < min2)
	    min2 = (check - solu.sol(0).pos(j)).abs2();
	}
	
	if(min2>minsol2){
	  // yes! we keep this solvent 
	  for(int k=0; k< num_atoms_per_solvent; k++){
	    solu.sol(0).addPos(solv.sol(0).pos(num_atoms_per_solvent * i + k));
          }
	}
      }
    }


    ostringstream title;
    title << "Solvating " << args["pos"];
    if (args.count("gather") != -1)
      title << " (gathered) ";
    
    title << " in " << args["solvent"] 
	  << endl;

    title << "Box dimensions (";
    if (boundary == truncoct) title << "truncated octahedron";
    else if (boundary == cubic) title << "cubic";
    else if (boundary == triclinic) title << "triclinic";
    else if (boundary == rectangular) title << "rectangular";
    else throw gromos::Exception("sim_box", "wrong boundary!!!");
    
    title << ") were ";
    if(!boxsize){
      title << "calculated from maximum" << endl;
      title << "solute atom-atom distance";
      if(args.count("rotate") == -1)
	title << " (not rotated):" << endl;
      else{
	if (boundary == rectangular)
	  title << "s (x, y, z) (after rotation):" << endl;
	else
	  title << " (after rotation):" << endl;
      }
      
      for(unsigned int i=0; i<max_dist.size(); ++i){
	int index = 2*(max_dist.size() - i - 1);
	title << "\t" << max_dist[i] << " between atoms "  
	      << as.mol(index)+1 << ":" << as.atom(index)+1 << " and " 
	      << as.mol(index+1)+1 << ":" << as.atom(index+1)+1 << endl;
      }
    }
    else 
      title <<"specified by user" << endl;
    if(numSolventAtoms){
      title << "System contained " << numSolventAtoms/num_atoms_per_solvent
	    << " solvent molecules" << endl;
    }
    title << "Added " << (solu.sol(0).numPos()-numSolventAtoms)
      / num_atoms_per_solvent
	  << " solvent molecules";

    // before writing out the data see if there were solute velocities read in
    // and put the missing solute velocities to zero if so
    bool writeVel = false;
    for(int mol = 0; mol < solu.numMolecules(); mol++) {
      if(solu.mol(mol).numVel() > 0) {
        writeVel = true;
      }
    }
    for(int solv = 0; solv < solu.numSolvents(); solv++) {
      if(solu.sol(solv).numVel() > 0) {
        writeVel = true;
      }
    }
    if (writeVel) {
      for (int mol = 0; mol < solu.numMolecules(); mol++) {
        if (solu.mol(mol).numVel() != solu.mol(mol).numPos()) {
          solu.mol(mol).initVel();
        }
      }
      for (int solv = 0; solv < solu.numSolvents(); solv++) {
        if (solu.sol(solv).numVel() < solu.sol(solv).numPos()) {
          int diff =solu.sol(solv).numPos() - solu.sol(solv).numVel();
          for (int i = 0; i < diff; i++) {
            solu.sol(solv).addVel(Vec(0.0, 0.0, 0.0));
          }
        }
      }
    }

    OutG96S oc(cout);      
    oc.select("ALL");
    solu.box().boxformat() = gcore::Box::genbox;
    
    oc.writeTitle(title.str());
    oc << solu;
    oc.close();
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
    
  }
  return 0;
}

void rotate_solute(System &sys, vector<double> &max, AtomSpecifier &as)
{
  max.resize(3);
  
  // this will be a two stage function.
  // First calculate the maximum distance between any two solute atoms
  max[0] =calc_max_size(sys, as, 3);
  
  // rotate the solute such that the atoms in as are along the z-axis
  Vec v=sys.mol(as.mol(0)).pos(as.atom(0)) 
    - sys.mol(as.mol(1)).pos(as.atom(1));
  double r = v.abs();
  double r_yz = sqrt(v[1]*v[1] + v[2]*v[2]);  
  
  // the rotation matrix is the product of two rotations
  // 1. around the x-axis by theta
  //    with sin(theta) = v[1]/r_yz; cos(theta) = v[2]/r_yz
  // 2. around the y-axis by phi
  //    with sin(phi) = v[0]/r; cos(phi) = r_yz / r
  if(r==0.0 || r_yz==0.0){
    throw gromos::Exception("sim_box",
			   "rotation of solute failed (z-dimension). Do you really want to use @rotate?");    
  }
  
  Matrix rot1(Vec( r_yz / r         ,  0         , v[0]/r ),
	      Vec(-v[0]*v[1]/r/r_yz ,  v[2]/r_yz , v[1]/r ),
	      Vec(-v[0]*v[2]/r/r_yz , -v[1]/r_yz , v[2]/r ));
  
  for(int m=0; m<sys.numMolecules(); m++)
    for(int a=0; a<sys.mol(m).numAtoms(); a++)
      sys.mol(m).pos(a) = rot1*sys.mol(m).pos(a);
  // and take along any solvent
  for(int i=0; i<sys.sol(0).numPos(); i++){
    sys.sol(0).pos(i) = rot1*sys.sol(0).pos(i);
  }
  
  // calculate the maximum distance in the x-y-plane
  max[1] =calc_max_size(sys, as, 2);

  // rotate the solute around the z-axis, such that the atoms in as are
  // along the y-axis, this is done by a rotation around psi with
  // sin(psi) = x/r_xy; cos(psi) = y/r_xy;
  v = sys.mol(as.mol(2)).pos(as.atom(2)) - sys.mol(as.mol(3)).pos(as.atom(3));
  double r_xy= sqrt(v[0]*v[0] + v[1] * v[1]);
  
  if(r_xy == 0.0){
    throw gromos::Exception("sim_box",
			    "rotation of solute failed (y-dimension). Do you really want to use @rotate?");
  }
  
  Matrix rot2(Vec( +v[1]/r_xy ,  v[0]/r_xy , 0),
	      Vec( -v[0]/r_xy ,  v[1]/r_xy , 0),
	      Vec( 0         ,  0         , 1));
  
  for(int m=0; m<sys.numMolecules(); m++)
    for(int a=0; a<sys.mol(m).numAtoms(); a++)
      sys.mol(m).pos(a) = rot2*sys.mol(m).pos(a);

  // and take along any solvent
  for(int i=0; i<sys.sol(0).numPos(); i++){
    sys.sol(0).pos(i) = rot2*sys.sol(0).pos(i);
  }
  
  // finally we calculate the maximum distance in the x-direction
  max[2] =calc_max_size(sys, as, 1);
}


double calc_max_size(System &sys, AtomSpecifier &as, int dim)
{
  // calculate the longest distance between solute atoms considering
  // the first dim dimensions.
  double max2=0.0;
  double d2=0;
  int max_m1=0, max_m2=0, max_a1=0, max_a2=0;
  
  for(int m1=0; m1 < sys.numMolecules(); m1++){
    for(int a1=0; a1 < sys.mol(m1).numAtoms(); a1++){
      Vec current=sys.mol(m1).pos(a1);
      for(int m2=m1; m2 < sys.numMolecules(); m2++){
	int start=0;
	if(m1==m2) start=a1;
	for(int a2=start; a2 < sys.mol(m2).numAtoms(); a2++){
	  d2=0.0;
	  for(int i=0; i<dim; i++){
	    d2 += (current[i] - sys.mol(m2).pos(a2)[i]) *
	      (current[i] - sys.mol(m2).pos(a2)[i]);
	  }
	  if(d2>max2){
	    max2=d2;
	    max_m1=m1;
	    max_m2=m2;
	    max_a1=a1;
	    max_a2=a2;
	  }
	}
      }
    }
  }
  as.addAtomStrict(max_m1, max_a1);
  as.addAtomStrict(max_m2, max_a2);
  
  return sqrt(max2);
}

