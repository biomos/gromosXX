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
 * @file bin_box.cc
 * Create a condensed phase system of two components
 */

/**
 * @page programs Program Documentation
 *
 * @anchor bin_box
 * @section bin_box Create a condensed phase system of two components
 * @author @ref dt
 * @date 7-6-07
 *
 * When simulating a molecular liquid, a starting configuration for the solvent
 * molecules has to be generated. To generate a starting configuration for the 
 * simulation of a binary mixture, the program bin_box can be used. A cubic box
 * is filled with solvent molecules by randomly distributing them on an evenly
 * spaced grid such that the total density of the box and the mole fractions of
 * the solvent components match the specified values. Note that as an
 * alternative, program @ref ran_box (see section V-2.10) can be used, which 
 * generates a starting configuration for the simulation of mixtures consisting
 * of an unlimited number of components, in which the molecules are oriented
 * randomly.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo1</td><td>&lt;molecular topology file (type 1)&gt; </td></tr>
 * <tr><td> \@pos1</td><td>&lt;coordinates for a single molecule of type 1&gt; </td></tr>
 * <tr><td> \@topo2</td><td>&lt;molecular topology file (type 2)&gt; </td></tr>
 * <tr><td> \@pos2</td><td>&lt;coordinates for a single molecule of type 2&gt; </td></tr>
 * <tr><td> \@nsm</td><td>&lt;number of molecules per dimension&gt; </td></tr>
 * <tr><td> \@densit</td><td>&lt;density of liquid (kg/m^3)&gt; </td></tr>
 * <tr><td> \@fraction</td><td>&lt;mole fraction of mixture component 1&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  bin_box
    @topo1    urea.top
    @pos1     urea.g96
    @topo2    h2o.top
    @pos2     h2o.g96
    @nsm      10
    @densit   888
    @fraction 0.2
 @endverbatim
 *
 * <hr>
 */
// buildbox.cc

#include <cassert>
#include <cstdlib>
#include <string>
#include <set>
#include <sstream>
#include <cmath>
#include <iostream>

#include "../src/args/Arguments.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Vec.h"
#include "../src/gromos/Exception.h"

using namespace gcore;
using namespace gio;
using namespace fit;
using namespace gmath;
using namespace args;
using namespace std;

class point
{
public:
  point(int i, int j, int k)
  {
    xi=i;
    yi=j;
    zi=k;
  }
  int xi;
  int yi;
  int zi;
  bool operator<(point const & p)const
  {
    if(xi<p.xi) return true;
    else if(xi==p.xi && yi < p.yi) return true;
    else if(xi==p.xi && yi==p.yi && zi < p.zi) return true;
    return false;
  }
};


int main(int argc, char **argv){

  Argument_List knowns;
  knowns << "topo1" << "pos1" << "topo2" << "pos2" << "nsm" << "densit"
         << "fraction";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo1    <molecular topology file (type 1)>\n";
  usage += "\t@pos1     <coordinates for a single molecule of type 1>\n";
  usage += "\t@topo2    <molecular topology file (type 2)>\n";
  usage += "\t@pos2     <coordinates for a single molecule of type 2>\n";
  usage += "\t@nsm      <number of molecules per dimension>\n";
  usage += "\t@densit   <density of liquid (kg/m^3)>\n";
  usage += "\t@fraction <mole fraction of mixture component 1>\n";

  try{
    Arguments args(argc, argv, knowns, usage);
    // set some values
    args.check("nsm",1);
    Arguments::const_iterator iter=args.lower_bound("nsm");
    int nsm3=atoi(iter->second.c_str());
    int nsm=nsm3*nsm3*nsm3;

    double densit=atof(args["densit"].c_str());
    double fraction=atof(args["fraction"].c_str());
    int nsm1=int(fraction*nsm);
    int nsm2=nsm-nsm1;
    
    // read topology
    args.check("topo1",1);
    InTopology it1(args["topo1"]);
    System smol1(it1.system());

    
    // read topology
    args.check("topo2",1);
    InTopology it2(args["topo2"]);
    System smol2(it2.system());

    // and calculate some more values
    double weight1=0, weight2=0;
    for(int i=0; i<smol1.numMolecules();i++)
      for(int j=0; j< smol1.mol(i).numAtoms();j++)
        weight1+=smol1.mol(i).topology().atom(j).mass();
    for(int i=0; i<smol2.numMolecules();i++)
      for(int j=0; j< smol2.mol(i).numAtoms(); j++)
	weight2+=smol2.mol(i).topology().atom(j).mass();
    double weight=nsm1*weight1 + nsm2*weight2;
    
    double vtot=(weight*1.66056)/densit;
    double box=pow(vtot,1.0/3.0);
    double box3=box/nsm3;
    Vec box32(box3/2.0, box3/2.0, box3/2.0);
    //cout << "total number of molecules " << nsm << endl;
    //cout << "number of type 1 " << nsm1 << endl;
    //cout << "number of type 2 " << nsm2 << endl;
    
    // read singe atom coordinates...
    InG96 ic;
    args.check("pos1",1);
    ic.open(args["pos1"]);
    ic >> smol1;
    ic.close();
    ic.open(args["pos2"]);
    ic >> smol2;
    ic.close();
        
    Vec rc1=PositionUtils::com(smol1)-box32;
    Vec rc2=PositionUtils::com(smol2)-box32;
    
    PositionUtils::translate(&smol1, -rc1);
    PositionUtils::translate(&smol2, -rc2);

     // new system
    System sys;
    //for(int i=0;i<3;i++){
    //  double *tmp = (double *) &sys.box()[i];
    //  *tmp = box;
    //}
    //initialize random seed
    srand(int(1000*densit));
   
    // set up a grid
    set<point> grid;
    for(int i=0; i<nsm3; i++)
      for(int j=0; j<nsm3; j++)
	for(int k=0; k<nsm3; k++)
	  grid.insert(point(i,j,k));
    
    //first, we do molecule 1
    for(int i=0; i< nsm1; i++){
      //get an random number
      int r=rand();
      long int itry=(r*grid.size())/RAND_MAX;
      set<point>::iterator iter=grid.begin();
      for(int i=0; i<itry; i++) iter++;
      if(iter==grid.end())
	throw gromos::Exception("bin_box", "unlikely but true: "
				"not enough grid points for random number");
    
      //cout << "itry " << itry << "\t" << ipos << "\t" << ix << "\t" << iy << "\t" << iz << endl;

      Vec shift(iter->xi*box3, iter->yi*box3, iter->zi*box3);
      PositionUtils::translate(&smol1, shift);
      for(int k=0; k< smol1.numMolecules(); k++)
	sys.addMolecule(smol1.mol(k));
      PositionUtils::translate(&smol1, -shift);
      grid.erase(iter);
    }
    //then, we do molecule 2
    if(int(grid.size())!=nsm2)
      throw gromos::Exception("bin_box", "the number of grid points left "
			      "after adding the first species is not the "
			      "same as required for the second species");
    // just fill in the remaining grid points
    for (set<point>::iterator iter = grid.begin(), to = grid.end();
            iter != to; ++iter) {

      Vec shift(iter->xi*box3, iter->yi*box3, iter->zi * box3);
      PositionUtils::translate(&smol2, shift);
      for (int k = 0; k < smol2.numMolecules(); k++)
        sys.addMolecule(smol2.mol(k));
      PositionUtils::translate(&smol2, -shift);
    }
    sys.box() = Box(gcore::Box::rectangular,
            double(box), double(box), double(box),
            90.0, 90.0, 90.0,
            0.0, 0.0, 0.0);


    // Print the new set to cout
    OutG96S oc;
    ostringstream os;
    os << "Bin_box generated a binary mixture of" << endl;
    os <<  nsm1 << " x "<<args["pos1"] << endl;
    os <<  nsm2 << " x "<<args["pos2"] << endl;
    os << "Density : " << densit << " kg/m^3";
    
    oc.open(cout);
    oc.writeTitle(string(os.str()));
    oc << sys;
    
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}




