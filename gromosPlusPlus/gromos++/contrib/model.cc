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
 * @file model.cc
 * models atom positions based on expressions
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor model
 * @section model models atom positions based on expressions
 * @author @ref ns
 * @date 11-12-2007
 *
 * Program models sets the atom position of given atoms to 
 * the result of an expression. The expression has to evalulate
 * to a vector value. See @ref utils::ExpressionProperty for 
 * details.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@pos</td><td>&lt;coordinate file&gt; </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms" to change&gt; </td></tr>
 * <tr><td> \@exprs</td><td>&lt;expressions for evert atoms&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  modeler
    @topo             ex.top
    @pbc              r
    @pos              ex.g96
    @atoms            1:1
    @exprs            expr%atom(1:1)+0.5*atom(1:2-3)/abs(atom(1:2-3))
 @endverbatim
 *
 * <hr>
 */


#include <cassert>
#include <cmath>
#include <cstdlib>
#include <map>
#include <map>
#include <map>
#include <map>
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
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/Value.h"
#include "../src/utils/Property.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/gromos/Exception.h"

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace gmath;
using namespace utils;
using namespace std;


int main(int argc, char **argv){
  Argument_List knowns;
  knowns << "topo" << "pbc" << "pos" << "atoms" << "exprs";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo           <molecular topology file>\n";
  usage += "\t@pbc            <boundary type> [<gathermethod>]\n";
  usage += "\t@pos            <coordinate file>\n";
  usage += "\t@atoms          <atoms to change>\n";
  usage += "\t@exprs          <expressions>\n";
  
 
  try{
    Arguments args(argc, argv, knowns, usage);
    
    InTopology it(args["topo"]);
    System sys(it.system());

    // open the file, get the system
    InG96 ic;
    ic.open(args["pos"]);
    ic.select("ALL");
    ic >> sys;
    ic.close();

    System refSys(it.system());

    Boundary *pbc = BoundaryParser::boundary(sys, args);
    Boundary::MemPtr gathmethod = GatherParser::parse(sys,refSys,args);
    // gather the system!
    (*pbc.*gathmethod)();

    
    AtomSpecifier atoms = AtomSpecifier(sys);
    {
      Arguments::const_iterator to = args.upper_bound("atoms");
      for(Arguments::const_iterator iter = args.lower_bound("atoms"); iter!=to;iter++)
        atoms.addSpecifier(iter->second);
    }  
    if(atoms.empty())
      throw gromos::Exception(argv[0], "No atoms specified.");
    
    PropertyContainer exprs(sys, pbc);
    {
      Arguments::const_iterator to = args.upper_bound("exprs");
      for(Arguments::const_iterator iter = args.lower_bound("exprs"); iter!=to;iter++)
        exprs.addSpecifier(iter->second);

      // check whether there are only expressions
      for(PropertyContainer::const_iterator it = exprs.begin(), to = exprs.end(); 
          it != to; ++it) {
        if ((*it)->type() != "Expression")
          throw gromos::Exception(argv[0], "@exprs does contain a non \"expr\" property.");
      }
      if (atoms.size() != int(exprs.size()))
        throw gromos::Exception(argv[0], "Please give an expression for every atom.");
    }

    // ready - model the atoms
    exprs.calc();
    
    ostringstream title;
    title << "model created atom positions for atoms:" << endl;

    // set the atom positions
    for(int i = 0; i < atoms.size(); ++i) {
      const Value posval = exprs[i]->getValue();
      
      Vec pos;
      try { // check for vector result
        atoms.pos(i) = posval.vec();
      } catch(...) {
        ostringstream msg;
        msg << "expression for atom " << i + 1 << " (" << atoms.toString(i) << ")"
            << " didn't evaluate to a vector.";
        throw gromos::Exception(argv[0], msg.str());
      }

      title << " - " << atoms.toString(i) << endl;
    }

    // write out the coordinates
    OutG96S oc(cout);
    oc.select("ALL");
    oc.writeTitle(title.str());
    oc << sys;
    oc.close();
  } catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
