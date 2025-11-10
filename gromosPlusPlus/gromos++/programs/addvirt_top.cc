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
 * @file addvirt_top.cc
 * create virtual atoms in a topology
 */

/**
 * @page programs Program Documentation
 *
 * @anchor addvirt_top
 * @section addvirt_top Add virtual atoms to a topology
 * @author @ref co
 * @date 29-10-2019
 *
 *
 * <hr>
 */


#include <cassert>
#include <cstdlib>
#include <set>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/gcore/System.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/OutTopology.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/VirtualAtom.h"
#include "../src/gcore/VirtualAtomType.h"
#include "../src/utils/Neighbours.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;


int main(int argc, char **argv){

  Argument_List knowns; 
  knowns << "topo" << "atomspec" << "param" << "hydrogens" << "addtype" << "addexclusions";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo  <molecular topology files>\n";
  usage += "\t[@atomspec <virtual atom specified>]\n";
  usage += "\t[@param    <iac> <charge> (for every virtual atom)]\n";
  usage += "\t[@hydrogens  <iac of united atoms to add H's for>]\n";
  usage += "\t[@addtype    <type> <dis1> <dis2>]\n";
  usage += "\t[@addexclusions]\n";

  try{
    Arguments args(argc, argv, knowns, usage);
    
    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());    
    GromosForceField gff(it.forceField());

    // first handle the virtual atom types
    for(Arguments::const_iterator iter=args.lower_bound("addtype"), 
          to=args.upper_bound("addtype"); iter!=to;++iter) {
      // read three numbers
      int type;
      double d1, d2;
      stringstream ss(iter->second);
      ss >> type;
      ++iter;
      if(iter==to){
	throw gromos::Exception("addvirt_top", 
            "@addtype needs to be followed by a multiple of 3 input parameters");
      }      
      ss.clear();
      ss.str(iter->second);
      ss >> d1; 
      ++iter;
      if(iter==to){
	throw gromos::Exception("addvirt_top", 
            "@addtype needs to be followed by a multiple of 3 input parameters");
      }      
      ss.clear();
      ss.str(iter->second);
      ss >> d2; 

      // check if exists
      if(gff.findVirtualAtomType(type)){
        cerr << "WARNING: Virtual atom type " << type 
             << " was already specified in the topology.\n"
             << "         Distances overwritten by user-specified input." << endl;
      }
      gff.addVirtualAtomType(VirtualAtomType(type, d1, d2));
    }
    gcore::VirtualAtoms vas = sys.vas();
    int numOrigVas = vas.numVirtualAtoms();

    // read atomspecifier
    utils::AtomSpecifier as(sys);
    for(Arguments::const_iterator iter=args.lower_bound("atomspec"), 
          to=args.upper_bound("atomspec"); iter!=to;++iter) {

      as.addSpecifier(iter->second);
    } 
    gcore::Exclusion e, e14;
    vas.addVirtualAtom(as, gff, -1, 0.0, e, e14);

    // give them parameters
    int count=numOrigVas;
    for(Arguments::const_iterator iter=args.lower_bound("param"),
        to=args.upper_bound("param"); iter!=to; ++iter) {
      istringstream is(iter->second);
      int iac = 1;
      if(!(is >> iac))
	throw gromos::Exception("addvirt_top", 
				"failed to read an integer for an IAC from @param "
				+ iter->second);
      iac--;
      ++iter;
      if(iter==args.upper_bound("param"))
        throw gromos::Exception("addvirt_top",
                                "specify both IAC and charge with @param");
      double charge=0;
      is.clear();
      is.str(iter->second);
      if(!(is >> charge))
 	throw gromos::Exception("addvirt_top", 
				"failed to read a double for a charge from @param "
				+ iter->second);
      if(count<vas.numVirtualAtoms()){
        vas.setIac(count, iac);
        vas.setCharge(count,charge);
      }
      else {
        cerr << "WARNING: Count mismatch between specified parameters and number of virtual atoms in @atomspec\n";
	cerr << "         Only used the first " << vas.numVirtualAtoms() - numOrigVas << " sets of parameters" << endl;
      }
      count++;
    } 

    if(count <vas.numVirtualAtoms()){
      cerr << "WARNING: Count mismatch between specified parameters and number of virtual atoms in @atomspec\n";
      cerr << "         Only set parameters for first " << count - numOrigVas << " virtual atoms" << endl;
    }
    // any explicit hydrogens
    std::set<int> uas;
    for(Arguments::const_iterator iter=args.lower_bound("hydrogens"), 
          to=args.upper_bound("hydrogens"); iter!=to;++iter) {

      int iac;
      istringstream is(iter->second);
      if(!(is >> iac))
	throw gromos::Exception("addvirt_top", 
				"failed to read an integer from @hydrogens"
				+ iter->second);

      uas.insert(iac -1);
    }
    // loop over all atoms in the system
    int totNumAtoms = 0;
    for(int m = 0; m < sys.numMolecules(); m++){
      for(int a = 0; a < sys.mol(m).numAtoms(); a++){
        if (uas.count(sys.mol(m).topology().atom(a).iac())) {
          utils::Neighbours n(sys, m, a);
          std::vector<int> n1, n2;
          gcore::Exclusion e, e14;
          switch(n.size()){
            case 1:
              {
              utils::Neighbours nn(sys, m, n[0]);
              if(nn.size() == 1){
                throw gromos::Exception("addvirt_top","cannot add hydrogens for a diatomic");
              }
              int k = nn[0];
              if(k==a) k=nn[1];
              n1.push_back(totNumAtoms + a);
              n1.push_back(totNumAtoms + n[0]);
              n1.push_back(totNumAtoms + k);
              vas.addVirtualAtom(sys, n1, 51, 
                                 gff.virtualAtomType(51).dis1(), 
                                 gff.virtualAtomType(51).dis2(),
                                 -1, 0.0, e, e14);
              vas.addVirtualAtom(sys, n1, 52, 
                                 gff.virtualAtomType(52).dis1(), 
                                 gff.virtualAtomType(52).dis2(),
                                 -1, 0.0, e, e14);
              vas.addVirtualAtom(sys, n1, 53, 
                                 gff.virtualAtomType(53).dis1(), 
                                 gff.virtualAtomType(53).dis2(),
                                 -1, 0.0, e, e14);
              }
              break;
            case 2:
              n1.push_back(totNumAtoms + a);
              n1.push_back(totNumAtoms + n[0]);
              n1.push_back(totNumAtoms + n[1]);
              vas.addVirtualAtom(sys, n1, 4, 
                                 gff.virtualAtomType(4).dis1(), 
                                 gff.virtualAtomType(4).dis2(),
                                 -1, 0.0, e, e14);
              n2.push_back(totNumAtoms + a);
              n2.push_back(totNumAtoms + n[1]);
              n2.push_back(totNumAtoms + n[0]);
              vas.addVirtualAtom(sys, n2, 4,
                                 gff.virtualAtomType(4).dis1(), 
                                 gff.virtualAtomType(4).dis2(),
                                 -1, 0.0, e, e14);
              break;
            case 3:
              n1.push_back(totNumAtoms + a);
              n1.push_back(totNumAtoms + n[0]);
              n1.push_back(totNumAtoms + n[1]);
              n1.push_back(totNumAtoms + n[2]);
              vas.addVirtualAtom(sys, n1, 1, 
                                 gff.virtualAtomType(1).dis1(), 
                                 gff.virtualAtomType(1).dis2(),
                                 -1, 0.0, e, e14);
              break;
            default:
              ostringstream error;
              error << "don't know how to add a virtual atom to atom " 
                    << a+1 << " in molecule " << m+1 
                    << " (option @hydrogens)";
              throw gromos::Exception("addvirt_top", error.str());
          }
        } 
      }
      totNumAtoms += sys.mol(m).numAtoms();
    }

    // add exclusions
    if(args.count("addexclusions") >= 0) {
      int numAtoms = 0; 
      for(int m = 0; m < sys.numMolecules(); m++){
        numAtoms+=sys.mol(m).numAtoms();
      }
      for(unsigned int i = numOrigVas; i < vas.numVirtualAtoms(); i++){
        // we add exclusions for the atoms that compose the virtual atoms
        for(unsigned int j=0; j < vas.atom(i).conf().size(); j++){
          sys.mol(vas.atom(i).conf().mol(j)).topology().atom(vas.atom(i).conf().atom(j)).exclusion().insert(i+numAtoms);
        } 
      } 
    }

    // put the VirtualAtoms back in the system
    sys.addVirtualAtoms(vas);

    // just to be sure, check the all the types for all the virtual atoms are
    // defined
    for(unsigned int i = 0; i< sys.vas().numVirtualAtoms(); i++){
      if(!gff.findVirtualAtomType(sys.vas().atom(i).type())){
        stringstream ss;
        ss << "Virtual atom with type " << sys.vas().atom(i).type() 
           << " defined, but the distances for this type are not included in "
           << "VIRTUALATOMTYPE block";
	throw gromos::Exception("addvirt_top",ss.str()); 
      }
    }
    OutTopology ot(cout);
  
    ostringstream title;
    title << it.title()
          << "addvirt_top added " << vas.numVirtualAtoms() - numOrigVas << " virtual atoms"; 
    ot.setTitle(title.str());
    ot.write(sys,gff);
  } catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

