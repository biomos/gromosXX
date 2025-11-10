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
 * @file top_to_build.cc
 * convert a topology to a building block
 */

/**
 * @page contrib Contrib Documentation
 *
 * @anchor top_to_build
 * @section top_to_build Convert a topology to a building block
 * @author @ref ns
 * @date 12.11.2010
 *
 * This program is used to convert a topology to a building block. In order to
 * carry out the mass mapping it reads the interaction function parameter file.
 * The molecule which is to be converted has to be specified.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@param</td><td>&lt;interaction function parameter file&gt; </td></tr>
 * <tr><td> \@mol</td><td>&lt;the molecule to convert&gt; </td></tr>
 * </table>
 *
 *
 * Example 1:
 * @verbatim
  add_atom
    @topo    ex.top
    @param   53A6.ifp
    @mol     1
 @endverbatim
 *
 * <hr>
 */
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>

#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InParameter.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/BuildingBlock.h"
#include "../src/gcore/BbSolute.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gromos/Exception.h"
#include "../src/gio/OutBuildingBlock.h"

using namespace gcore;
using namespace gio;
using namespace args;
using namespace std;


int main(int argc, char *argv[]){

  Argument_List knowns;
  knowns << "topo" << "param" << "mol";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@param <interaction function parameter file>\n";
  usage += "\t@mol   <the molecule to extract>\n";
  
  try{
    Arguments args(argc, argv, knowns, usage);
    
    InTopology it(args["topo"]);
    InParameter ip(args["param"]);
    int mol = args.getValue<int>("mol") - 1;

    System sys(it.system());
    if(mol < 0 || mol >= sys.numMolecules()){
      throw gromos::Exception(argv[0], "molecule out of range.");
    }

    
    MoleculeTopology nbb=sys.mol(mol).topology();
    BbSolute bb(nbb);
    for(int i = 0; i < nbb.numAtoms(); ++i) {
      int massCode = ip.forceField().findMassType(nbb.atom(i).mass());
      if (massCode < 0) {
        ostringstream msg;
        msg << "The mass " << nbb.atom(i).mass() << " has no corresponding "
                << "mass type code in the parameter file. Add it and repeat.";
        throw gromos::Exception(argv[0], msg.str());
      }
      bb.atom(i).setMass(massCode);
    }

    OutBuildingBlock obb(cout);
    obb.writeSingle(bb, OutBuildingBlock::BBTypeSolute);
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}





