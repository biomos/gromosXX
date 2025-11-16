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
 * @file copy_box.cc
 * Repeats a simulation box in a given direction
 */

/**
 * @page programs Program Documentation
 *
 * @anchor copy_box
 * @section copy_box Repeats a simulation box in a given direction
 * @author @ref co
 * @date 8-6-07
 *
 * Program copy box can be used to duplicate the size of a system in the x, y 
 * or z direction (or k,l,m for triclinic boxes). This is especially convenient
 * if one wants to double the size of a system under periodic boundary
 * conditions in which the central box has a rectangular or triclinic shape. If
 * one wants perform more elaborate transformations, the program @ref cry might
 * be of use (see section V-2.17). Note that program @ref com_top (see section
 * V-2.2) can be useful to additionally duplicate the solute block in the 
 * topology.
 * Note that the \@pbc flag is optional. Only if this flag is given, gathering
 * of the molecules will be performed before copying the box!
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pos</td><td>&lt;input coordinate file&gt; </td></tr>
 * <tr><td> \@dir</td><td>&lt;coordinate to duplicate: x/y/z/k/l/m&gt; </td></tr>
 * <tr><td> [\@pbc</td><td>&lt;boundary type&gt; </td></tr>;] 
 * </table>
 *
 *
 * Example:
 * @verbatim
  copy_box
    @topo ex.top
    @pos  exref.coo
    @dir  x
 @endverbatim
 *
 * <hr>
 */


#include <cassert>
#include <cstdlib>
#include <string>
#include <sstream>
#include <cmath>
#include <iostream>

#include "../src/args/Arguments.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Vec.h"
#include "../src/bound/Boundary.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace fit;
using namespace gmath;
using namespace bound;
using namespace args;
using namespace utils;

int main(int argc, char **argv){

  Argument_List knowns;
  knowns << "topo" << "pos" << "dir" << "pbc";  

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo <molecular topology file>\n";
  usage += "\t@pos  <input coordinate file>\n";
  usage += "\t@dir  <coordinate to duplicate: x/y/z/k/l/m>\n";
  usage += "\t[@pbc <boundary type> [<gathermethod>] ]\n";
  

  try{
    Arguments args(argc, argv, knowns, usage);

    ostringstream title;
    
    // set some values
    args.check("dir",1);
    Arguments::const_iterator iter=args.lower_bound("dir");
    char dir=iter->second.c_str()[0];
    
    // read topology
    args.check("topo",1);
    InTopology it(args["topo"]);
    System sys(it.system());

    System refSys(it.system());
    
    // read singe atom coordinates...
    InG96 ic;
    ic.open(args["pos"]);
    ic.select("ALL");
    ic >> sys;
    title << ic.title();    
    ic.close();


    // Check if @pbc is given
    if(args.count("pbc")>0){
      Boundary *pbc = BoundaryParser::boundary(sys,args);
      Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);
      (*pbc.*gathmethod)();
      if(dir=='x'||dir=='y'||dir=='z'){
        if(pbc->type()!='r'){
          std::cerr << "WARNING: in copy_box, you want to copy ";
          std::cerr << "a non-rectangular box along x, y, or z!";
          std::cerr << endl; 
        }
      }
      if(dir=='k'||dir=='l'||dir=='m'){
        if(pbc->type()!='c'){
          std::cerr << "WARNING: in copy_box, you want to copy ";
          std::cerr << "a non-triclinic box along k, l, or m!";
          std::cerr << endl;
        }
      }
    }
  
    //calculate shifting vector
    Vec shift;
    if(dir=='x') shift = Vec(sys.box().K()[0],0,0);
    else if(dir=='y') shift = Vec(0,sys.box().L()[1],0);
    else if(dir=='z') shift = Vec(0,0,sys.box().M()[2]);
    else if(dir=='k') shift = sys.box().K();
    else if(dir=='l') shift = sys.box().L();
    else if(dir=='m') shift = sys.box().M();
    else throw gromos::Exception("copy_box", 
         "invalid direction specified, select x,y or z (or k,l,m for triclinic boxes)");
    
    //copy and move all molecules
    System sy2(sys);
    PositionUtils::translate(&sy2, shift);
    for(int i=0; i< sy2.numMolecules();i++){
      sys.addMolecule(sy2.mol(i));
    }
    for(int i=0; i< sy2.sol(0).numPos();i++){
      sys.sol(0).addPos(sy2.sol(0).pos(i));
    }

    if(dir=='x')      sys.box().K()[0] += shift[0];
    else if(dir=='y') sys.box().L()[1] += shift[1];
    else if(dir=='z') sys.box().M()[2] += shift[2];
    else if(dir=='k') sys.box().K() += shift;
    else if(dir=='l') sys.box().L() += shift;
    else if(dir=='m') sys.box().M() += shift;
    // Print the new set to cout
    OutG96S oc;
    title << "\nCopy_box: " << args["pos"] << " duplicated in " 
          << dir << "-direction";
    
    oc.open(cout);
    oc.select("ALL");
    
    oc.writeTitle(title.str());
    oc << sys;
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}




