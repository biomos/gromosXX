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
 * @file mod_top.cc
 * Modify the order of atoms in a topology
 */

/**
 * @page programs Program Documentation
 *
 * @anchor mod_top
 * @section mod_top modify the order of atoms in a topology
 * @author @ref mp
 * @date 21-3-2017
 *
 * Program mod_top modifies the order in which atoms occur in the topology. This
 * can be useful as a first step before creating a perturbation topology using 
 * @ref make_pt_top from two already existing topologies.
 *
 * The user has to give a list of pairs of atom indices and the atom at the first
 * index position will be moved to the second index position. Multiple pairs can be given, 
 * but this requires extra care as the procedure is applied sequentially, i.e. 
 * the second move will be applied on an already modified topology.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@move</td><td>&lt;oldpos newpos (pairs of atom indices)&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  mod_top
    @topo   ex.top
    @move   7 2 8 15 
 @endverbatim
 *
 * <hr>
 */
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/OutTopology.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/LinearTopology.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;

int main(int argc, char *argv[]){
  Argument_List knowns;
  knowns << "topo" << "move";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo  <molecular topology file>\n";
  usage += "\t@move <oldpos newpos (pairs of atom indices)>\n";

  try{
    Arguments args(argc, argv, knowns, usage);

    InTopology it(args["topo"]);
    System sys(it.system());
    
    cerr << "# NOTE: You need to manually take care of any changes in the charge groups!\n\n";
    
    // Parse atom indices to be moved
    std::vector<std::pair<int,int> > pvec;
    for(Arguments::const_iterator iter=args.lower_bound("move"), 
          to=args.upper_bound("move"); iter!=to;++iter) {
        int a, b;
        a=atoi(iter->second.c_str());
        iter++;
        if (iter != to)   b=atoi(iter->second.c_str());
        else throw(gromos::Exception("mod_top", "@move : uneven number of values, but we need pairs!"));
        std::pair<int,int> p(a,b);
        pvec.push_back(p);
    }

    // create a linear topology
    gcore::LinearTopology lt(sys);
    
    // move atoms
    lt.moveAtoms(pvec);

    // and parse the linearized thing back into a topology (which might
    // have a zillion molecules, because bonds have been removed)
    System syo = lt.parse();
    
    // take the old solvent
    syo.addSolvent(sys.sol(0));
    
    // set the temperature and pressure groups
    {
      int a = 0;
      for (int m = 0; m < syo.numMolecules(); ++m) {
        a += syo.mol(m).numAtoms();
        syo.addPressureGroup(a);
        syo.addTemperatureGroup(a);
      }
      // compare the current number with the previous one and warn if it changed for other reasons than deleted molecuels
      if (sys.numMolecules() != sys.numPressureGroups()) {
        if (sys.numPressureGroups() != syo.numPressureGroups()) {
          cerr << "WARNING: The number of pressure groups has changed. manual check recommended.\n";
          cerr << "         Number of pressure groups: " << sys.numPressureGroups() << " (" << args["topo"] << ")" << endl
                  << "                                     " << syo.numPressureGroups() << " (reduced topology)\n";
        }
      }
      if (sys.numMolecules() != sys.numTemperatureGroups()) {
        if (sys.numTemperatureGroups() != syo.numTemperatureGroups()) {
          cerr << "WARNING: The number of temperature groups has changed. manual check recommended.\n";
          cerr << "         Number of temperature groups: " << sys.numTemperatureGroups() << " (" << args["topo"] << ")" << endl
                  << "                                        " << syo.numTemperatureGroups() << " (reduced topology)\n";
        }
      }
    }
    
    // and write out the new topology
    OutTopology ot(cout);
    ostringstream os;
    os << "Modified topology based on " << args["topo"] << endl;
    os << "changed atom positions: ";
    for(unsigned int i=0; i< pvec.size(); i++)
      os << pvec[i].first << "->" << pvec[i].second;
    os  << endl;
    ot.setTitle(os.str());
    
    ot.write(syo, it.forceField());
    
    return 0;
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}



