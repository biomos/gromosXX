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
 * @file amber2gromos.cc
 * Translates an AMBER topology into GROMOS format
 */

/**
 * @page programs Program Documentation
 *
 * @anchor amber2gromos
 * @section amber2gromos Translates an AMBER topology into GROMOS format
 * @author @ref sr
 * @date 4.3.16
 *
 * Description
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@ambertop</td><td>&lt;AMBER molecular topology file &gt; </td></tr>
 * <tr><td> \@solvent</td><td>&lt;GROMOS topology file with solvent&gt; </td></tr>
 * <tr><td> \[@ljscaling</td><td>&lt;scaling factor for LJ parameters (default: 2.0)]&gt; </td></tr>
 * <tr><td> \[@atomic_chargegroups</td><td>&lt;each atom is assigned to its own charge group (default: 0)]&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  amber2gromos
    @ambertop             ligand.prm
    @solvent              spc.top
    @ljscaling            2.0
    @atomic_chargegroups  0
    @chargegroups         <path to chargroup file>
 @endverbatim
 *
 * <hr>
 */

#include <cstdlib>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>

#include "../src/gromos/Exception.h"
#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/OutTopology.h"
#include "../src/gio/InAmberTopology.h"
#include "../src/gio/InAmberTopology.h"
#include "../src/gio/InChargeGroups.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/LinearTopology.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;

int main(int argc, char *argv[]){
   
  Argument_List knowns;
  knowns << "ambertop" << "solvent" << "ljscaling" <<"atomic_chargegroups" <<"chargegroups";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@ambertop  <AMBER molecular topology file>\n";
  usage += "\n\t@solvent   <GROMOS topology file with solvent>\n";
  usage += "\n\t[@ljscaling <scaling factor for LJ parameters (default: 2.0)>]\n";
  usage += "\n\t[@atomic_chargegroups <assign each atom to its own chargegroup (default: 0)>]\n";
  usage += "\n\t[@chargegroups <path to chargroup file>]\n";
  
  try{
    Arguments args(argc, argv, knowns, usage);

    // our AMBER topology         
    AmberTopology ambertop(args["ambertop"].c_str());
    double ljscaling = args.getValue<double>("ljscaling", false, 2.0);
    bool atomic_chargegroups = args.getValue<bool>("atomic_chargegroups", false, 0);

    // read in the AMBER topology 
    LinearTopology lt;
    ambertop.parseFile(lt, 1.0/ljscaling, atomic_chargegroups);

    //cerr << lt.atoms().size() << endl;
    
    // parse everything into a system  
    System sys;
    lt.parse(sys);

    // ChargeGroups
    // define new charge groups
    string chargeGroupsPath = args.getValue<string>("chargegroups", false, "none");
    
    if(chargeGroupsPath == "none"){ //
        string msg = "WARNING!\t\t ChargeGroupfile was not provided";
        cerr << msg << endl;
    } 
    else {
        InChargeGroups chargeGroups(chargeGroupsPath);
        sys = chargeGroups.mapChargeGroupsOnResidues(sys);
    }    

    // solvent
    // we need to have a solvent, although it will not be used
    // read the GROMOS topology
    //repair solvent Errors are triggering segmentation faults!
    string solventTopPath = args.getValue<string>("solvent", false, "none");
    
    if(solventTopPath == "none"){
        string msg = "No solvent topology was provided! (Please add @solvent argument with path to solvent.top\n";
        throw gromos::Exception("AMBER2Gromos",msg);
    }
    else{
        InTopology it(solventTopPath);
        sys.addSolvent(it.system().sol(0));
        //TODO: add Solvent params!
        //SolventTopology st(solventTopPath);
        //sys.addSolvent(Solvent(st));
    }
    
    // set the hydrogens
    for (int m = 0; m < sys.numMolecules(); ++m) {
      sys.mol(m).topology().clearH();
      sys.mol(m).topology().setHmass(1.008);
    }

    // set the temperature and pressure groups
    int totNumAt = 0;
    for (int i = 0; i < sys.numMolecules(); ++i) {
      totNumAt += sys.mol(i).numAtoms();
      sys.addTemperatureGroup(totNumAt);
      sys.addPressureGroup(totNumAt);
    }

    // write the topology
    OutTopology ot(cout);
    ostringstream title;
    title << " AMBER topology translated into GROMOS";
    
    ot.setTitle(title.str());
    ot.write(sys, ambertop.d_gff);
    cout.flush();
    exit(0);
  }
  catch(gromos::Exception e){
    cerr << "FOUND ERROR!: ";
    cerr << e.what() << endl;
    exit(1);
  }
}


