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
 * @file mut_top.cc
 * Mutate a specific residue in a topology
 */

/**
 * @page programs Program Documentation
 *
 * @anchor mut_top
 * @section mut_top Mutate a specific residue in a molecular topology
 * @author @ref co
 * @date 18-1-12
 *
 * For complex topologies, it may not be simple to regenerate a mutated version
 * of it. This program allows the user to replace one or more indicated residue
 * by a (single) different one in the molecular topology. 
 *
 * The program does not work for residues with special connectivities, like 
 * proline or the first or last residue in a chain. Also, if the residue is 
 * cross-linked to other parts of the topology (e.g. Cys-bridges), things may
 * not lead to the desired result. Always check if your output makes sense!
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@build</td><td>&lt;building block file&gt; </td></tr>
 * <tr><td> \@param</td><td>&lt;parameter file&gt; </td></tr>
 * <tr><td> \@mutate</td><td>&lt;atomspecifier for the residue to replace (all atoms)&gt; </td></tr>
 * <tr><td> \@to</td><td>&lt;name of the residue (building block to replace it with&gt; </td></tr></table>
 *
 *
 * Example:
 * @verbatim
  mut_top
    @topo   ex.top
    @build  54a8.mtb
    @param  54a8.ifp
    @mutate 1:res(4:a)
    @to     ARG
 @endverbatim
 *
 * <hr>
 */
#include <cassert>
#include <iostream>
#include <sstream>
#include <set>
#include <map>
#include <string>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InBuildingBlock.h"
#include "../src/gio/InParameter.h"
#include "../src/gio/OutTopology.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/CrossDihedral.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/LinearTopology.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/gcore/BuildingBlock.h"
#include "../src/gcore/BbSolute.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;

int main(int argc, char *argv[]){
  Argument_List knowns;
  knowns << "topo" << "build" << "param" << "mutate" << "to";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t@build  <building block file>\n";
  usage += "\t@param  <parameter file>\n";
  usage += "\t@mutate <atomspecifier for the residue to replace (all atoms)>\n";
  usage += "\t@to     <name of the residue (building block to replace it with>\n";

  try{
    Arguments args(argc, argv, knowns, usage);

    InTopology it(args["topo"]);
    System sys(it.system());

    // read in the building block file
    InBuildingBlock ibb(args["build"]);
    BuildingBlock mtb(ibb.building());

    // read in the parameter file (we need it to translate the masstypes)
    InParameter ip(args["param"]);
    GromosForceField gff(ip.forceField());
    
    // define the selected block
    int index=mtb.findBb(args["to"]);
    if(index==0)
      throw gromos::Exception("mut_top",
	      "Building block " + args["block"] + " not found" );
    BbSolute bb;
    
    if(index>0)
      bb=mtb.bb(index-1);
    else{
      throw gromos::Exception("mut_top",
	    "Building block " + args["block"] + " is an MTBUILDBLEND group, which is not supported" );
    }

    // create a linear topology
    gcore::LinearTopology lt(sys);
    
    // and an empty, new, one
    gcore::LinearTopology ltnew;
    
    // find the residue to mutate
    utils::AtomSpecifier as(sys);
    for(Arguments::const_iterator iter=args.lower_bound("mutate"), 
          to=args.upper_bound("mutate"); iter!=to;++iter) {
      as.addSpecifier(iter->second);
    }
    
    // do the checks for known things that won't work
    if(as.resnum(0)==0){
      throw gromos::Exception("mut_top",
			       "This program cannot be used to mutate the first residue");
    }
    if(as.resnum(as.size()-1)==sys.mol(as.mol(as.size()-1)).topology().numRes()-1){
      throw gromos::Exception("mut_top",
			       "This program cannot be used to mutate the last residue of a molecule");
    }
    // ugly hardcoded names
    if(args["to"]=="PRO" || args["to"]=="HYP"){
      throw gromos::Exception("mut_top",
			      "This program cannot be used to mutate to a residue that does not have the same connectivity in the backbone");
    }
    
    
    // do some math
    int unchanged = as.gromosAtom(0);
    int num_to_add = bb.numAtoms();
    int num_to_remove = as.size();
    int diff = num_to_add - num_to_remove;
    std::vector<int> ren;
    for(int i=0; i< unchanged+3; i++) 
      ren.push_back(i);
    for(int i=0; i< num_to_remove-5; i++) 
      ren.push_back(-1);
    for(unsigned int i=unchanged+num_to_remove-2; i < lt.atoms().size(); i++)
      ren.push_back(i+diff);

    // now we can start building up the topology
    // first the atoms until unchanged+2; (we leave the CA for the new one)
    for(int i=0; i< unchanged; i++){
      AtomTopology aa(lt.atoms()[i]);
      Exclusion e;
      for(int j=0; j < lt.atoms()[i].exclusion().size(); j++){
	if(ren[lt.atoms()[i].exclusion().atom(j)]!=-1)
	  e.insert(ren[lt.atoms()[i].exclusion().atom(j)]);
      }
      aa.setExclusion(e);
      ltnew.addAtom(aa);
      // the residue count is still the same here
      ltnew.setResNum(ltnew.atoms().size()-1, lt.resMap()[i]);
      ltnew.setResName(lt.resMap()[i], lt.resNames()[lt.resMap()[i]]);
      
    }
    //std::cout << "added " << ltnew.atoms().size() << std::endl;
    //std::cout << "num_to_add " << num_to_add << std::endl;

    // we overwrite the exclusions of the last numPexcl atoms, using the
    // previous exclusions in the latest building block.
    for(int i=0; i< bb.numPexcl(); i++){
      Exclusion e;
      for(int j=0; j< bb.pexcl(i).size(); j++){
	e.insert(bb.pexcl(i).atom(j) + unchanged);
      }
      ltnew.atoms()[ltnew.atoms().size()-bb.numPexcl()+i].setExclusion(e);
    }
		 
    // what is the highest residue number
    int rescount = ltnew.resMap()[ltnew.atoms().size()-1]+1;
    
    // add the new atoms, starting from CA
    for(int i=0; i< num_to_add-2; i++){
      //std::cout << i << std::endl;
      
      AtomTopology aa(bb.atom(i));
      aa.setMass(gff.findMass(int(bb.atom(i).mass())));
      // hardcoded check for H
      if(aa.mass()==1.008) aa.setH(true);
      
      Exclusion e;
      for(int j=0; j < bb.atom(i).exclusion().size(); j++){
	//std::cout << "exlcusion " << bb.atom(i).exclusion().atom(j) << "-> " <<  bb.atom(i).exclusion().atom(j)+unchanged << std::endl;
	
	e.insert(bb.atom(i).exclusion().atom(j)+unchanged);
	
      }
      aa.setExclusion(e);
      ltnew.addAtom(aa);
      ltnew.setResNum(ltnew.atoms().size()-1,  rescount);
      ltnew.setResName(rescount, args["to"]);
    }
    //std::cout << "added " << ltnew.atoms().size() << std::endl;


   // and the remaining ones
    for(unsigned int i=unchanged+num_to_remove-2; i< lt.atoms().size(); i++){
      AtomTopology aa(lt.atoms()[i]);
      Exclusion e;
      for(int j=0; j < lt.atoms()[i].exclusion().size(); j++){
	if(ren[lt.atoms()[i].exclusion().atom(j)]!=-1)
	  e.insert(ren[lt.atoms()[i].exclusion().atom(j)]);
      }
      aa.setExclusion(e);
      ltnew.addAtom(aa);
      if(lt.resMap()[i] != lt.resMap()[i-1]){
	rescount++;
	ltnew.setResName(rescount, lt.resNames()[lt.resMap()[i]]);
      }
      
      ltnew.setResNum(ltnew.atoms().size()-1,  rescount);
    }
    
    //std::cout << "ltnew has " << ltnew.atoms().size() << " atoms\n";
    
    // now we can add all bonds
    for(std::set<Bond>::iterator iter=lt.bonds().begin(), to=lt.bonds().end();
	iter!=to; iter++){
      if(ren[ (*iter)[0] ] != -1 && ren[ (*iter)[1] ] != -1){
	  
	Bond b(ren[(*iter)[0]],ren[(*iter)[1]]);
	b.setType(iter->type());
	ltnew.addBond(b);
      }
    }
    // and the new ones
    BondIterator bi(bb);
    for(;bi; ++bi){
      Bond b(bi()[0]+unchanged, bi()[1]+unchanged);
      b.setType(bi().type());
      ltnew.addBond(b);
    }
    //std::cout << "numBonds " << ltnew.bonds().size() << std::endl;
    
    // now we can add all angles
    for(std::set<Angle>::iterator iter=lt.angles().begin(), to=lt.angles().end();
	iter!=to; iter++){
      if(ren[(*iter)[0]] != -1 && 
	 ren[(*iter)[1]] != -1 && 
	 ren[(*iter)[2]] != -1 ){
	  
	Angle b(ren[(*iter)[0]],ren[(*iter)[1]], ren[(*iter)[2]]);
	b.setType(iter->type());
	ltnew.addAngle(b);
      }
    }
    // and the new ones
    AngleIterator ai(bb);
    for(;ai; ++ai){
      Angle b(ai()[0]+unchanged, ai()[1]+unchanged, ai()[2] +unchanged);
      b.setType(ai().type());
      ltnew.addAngle(b);
    }
    //std::cout << "numAngles " << ltnew.angles().size() << std::endl;
    
    // now we can add all Impropers
    for(std::set<Improper>::iterator iter=lt.impropers().begin(), to=lt.impropers().end();
	iter!=to; iter++){
      if(ren[(*iter)[0]] != -1 && 
	 ren[(*iter)[1]] != -1 && 
	 ren[(*iter)[2]] != -1 &&
	 ren[(*iter)[3]] != -1){
	  
	Improper b(ren[(*iter)[0]],ren[(*iter)[1]], ren[(*iter)[2]],ren[(*iter)[3]] );
	b.setType(iter->type());
	ltnew.addImproper(b);
      }
    }
    // and the new ones
    ImproperIterator ii(bb);
    for(;ii; ++ii){
      Improper b(ii()[0]+unchanged, ii()[1]+unchanged, ii()[2] +unchanged,ii()[3] +unchanged );
      b.setType(ii().type());
      ltnew.addImproper(b);
    }
    //std::cout << "numImpropers " << ltnew.impropers().size() << std::endl;
    
    // now we can add all Dihedrals
    for(std::set<Dihedral>::iterator iter=lt.dihedrals().begin(), to=lt.dihedrals().end();
	iter!=to; iter++){
      if(ren[(*iter)[0]] != -1 && 
	 ren[(*iter)[1]] != -1 && 
	 ren[(*iter)[2]] != -1 &&
	 ren[(*iter)[3]] != -1){
	  
	Dihedral b(ren[(*iter)[0]],ren[(*iter)[1]], ren[(*iter)[2]],ren[(*iter)[3]] );

	b.setType(iter->type());
	//std::cout << "adding dihedral " << b[0] << " " << b[1] << " " << b[2] << " " << b[3] << " : " << b.type() << std::endl;
	ltnew.addDihedral(b);
      }
      // special case: the first atom is still there this will allow us to remove more than one residue by a single new one
      if(ren[(*iter)[0]] == -1 &&
	 ren[(*iter)[1]] != -1 && 
	 ren[(*iter)[2]] != -1 &&
	 ren[(*iter)[3]] != -1){
	
	int corr0 = 0;
	for (std::set<gcore::Bond>::const_iterator itbond = ltnew.bonds().begin();
	     itbond != ltnew.bonds().end(); ++itbond)
	  if ((*itbond)[1] == ren[(*iter)[1]]) corr0 = (*itbond)[0];
      

	Dihedral b(corr0,ren[(*iter)[1]], ren[(*iter)[2]],ren[(*iter)[3]] );
	b.setType((*iter).type());
	ltnew.addDihedral(b);
	//std::cerr << "adding dihedral " << b[0] << " " << b[1] << " " << b[2] << " " << b[3] << " : " << b.type()  << std::endl;
      }
      
    }
    // and the new ones
    DihedralIterator di(bb);
    for(;di; ++di){
      // we skip the ones to atom -2, it is already included and otherwise we have to find the atom that it refers to
      if(di()[0] != -3 && di()[3] != -3){
	  
	Dihedral b(di()[0]+unchanged, di()[1]+unchanged, di()[2] +unchanged,di()[3] +unchanged );
	b.setType(di().type());
	ltnew.addDihedral(b);
	//std::cout << "adding dihedral " << di()[0] << " " << di()[1] << " " << di()[2] << " " << di()[3] << " : " << di().type()  << std::endl;
	
      }
      else if(di()[0]==-3){
	//cerr << "a peptide bond" << std::endl;
	
	int corr0 = 0;
	for (std::set<gcore::Bond>::const_iterator iter = ltnew.bonds().begin();
	     iter != ltnew.bonds().end(); ++iter)
	  if ((*iter)[1] == di()[1] + unchanged) corr0 = (*iter)[0] + 3;
      

	Dihedral b(di()[0] + corr0, di()[1] + unchanged, di()[2] + unchanged, di()[3] + unchanged);
	b.setType(di().type());
	ltnew.addDihedral(b);
	//std::cerr << "adding dihedral " << di()[0] << " " << di()[1] << " " << di()[2] << " " << di()[3] << " : " << di().type()  << std::endl;
	//std::cerr << "adding dihedral " << b[0] << " " << b[1] << " " << b[2] << " " << b[3] << " : " << b.type()  << std::endl;
      }
      else if(di()[3]==-3){
	int corr3 = 0;
	for (std::set<gcore::Bond>::const_iterator iter = ltnew.bonds().begin();
	     iter != ltnew.bonds().end(); ++iter)
	  if ((*iter)[1] == di()[2] + unchanged) corr3 = (*iter)[0] + 3;
	
	
	Dihedral b(di()[0] + unchanged, di()[1] + unchanged, di()[2] + unchanged, di()[3] + corr3);
	b.setType(di().type());
	ltnew.addDihedral(b);
      }
      
    }
    
    //std::cout << "numDihedrals " << ltnew.dihedrals().size() << std::endl;
    
    // now we can add all CrossDihedrals
    for(std::set<CrossDihedral>::iterator iter=lt.crossdihedrals().begin(), to=lt.crossdihedrals().end();
	iter!=to; iter++){
      if(ren[(*iter)[0]] != -1 && 
	 ren[(*iter)[1]] != -1 && 
	 ren[(*iter)[2]] != -1 &&
	 ren[(*iter)[3]] != -1 &&
	 ren[(*iter)[4]] != -1 && 
	 ren[(*iter)[5]] != -1 && 
	 ren[(*iter)[6]] != -1 &&
	 ren[(*iter)[7]] != -1){
	  
	CrossDihedral b(ren[(*iter)[0]],ren[(*iter)[1]], ren[(*iter)[2]],ren[(*iter)[3]],ren[(*iter)[4]],ren[(*iter)[5]], ren[(*iter)[6]],ren[(*iter)[7]] );
	b.setType(iter->type());
	ltnew.addCrossDihedral(b);
      }
    }
    // and the new ones
    CrossDihedralIterator ci(bb);
    for(;ci; ++ci){
      CrossDihedral b(ci()[0]+unchanged, ci()[1]+unchanged, ci()[2] +unchanged,ci()[3] +unchanged,ci()[4]+unchanged, ci()[5]+unchanged, ci()[6] +unchanged,ci()[7] +unchanged );
      b.setType(ci().type());
      ltnew.addCrossDihedral(b);
    }
    //std::cout << "numCrossDihedrals " << ltnew.crossdihedrals().size() << std::endl;
    
    // calculate the new 1,4 interactions
    ltnew.get14s();
    
    System syo = ltnew.parse(); 
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
          cerr << "         Number of pressure greoups: " << sys.numPressureGroups() << " (" << args["topo"] << ")" << endl
                  << "                                     " << syo.numPressureGroups() << " (reduced topology)\n";
        }
      }
      if (sys.numMolecules() != sys.numTemperatureGroups()) {
        if (sys.numTemperatureGroups() != syo.numTemperatureGroups()) {
          cerr << "WARNING: The number of temperature groups has changed. manual check recommended.\n";
          cerr << "         Number of temperature greoups: " << sys.numTemperatureGroups() << " (" << args["topo"] << ")" << endl
                  << "                                        " << syo.numTemperatureGroups() << " (reduced topology)\n";
        }
      }
    }
    
    // and write out the new topology
    vector<string> as_str = as.toString();

    OutTopology ot(cout);
    ostringstream os;
    os << "Mutated topology based on " << args["topo"] << endl;
    os << "in which atoms ";
    for(unsigned int i=0; i< as_str.size(); i++)
      os << as_str[i] << " ";
    os << "are replaced by residue " << args["to"] << endl;
    os << "WARNING: This program has not been thoroughly tested: check your output!" << endl;
    
    ot.setTitle(os.str());
    
    ot.write(syo, it.forceField());
   
    /*
    

    // Parse atom specifiers
    utils::AtomSpecifier as(sys);
    for(Arguments::const_iterator iter=args.lower_bound("atoms"), 
          to=args.upper_bound("atoms"); iter!=to;++iter) {
      as.addSpecifier(iter->second);
    }
    vector<string> as_str = as.toString();

    // flag all atoms that are not in the list with a negative iac
    for(int m=0; m<sys.numMolecules(); m++){
      for(int a=0; a<sys.mol(m).numAtoms(); a++){
	if(as.findAtom(m,a)==-1) sys.mol(m).topology().atom(a).setIac(-1);
      }
    }

    // create a linear topology
    gcore::LinearTopology lt(sys);

    // remove all flagged atoms
    lt.removeAtoms();
    
    // calculate the new 1,4 interactions
    lt.get14s();

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
          cerr << "         Number of pressure greoups: " << sys.numPressureGroups() << " (" << args["topo"] << ")" << endl
                  << "                                     " << syo.numPressureGroups() << " (reduced topology)\n";
        }
      }
      if (sys.numMolecules() != sys.numTemperatureGroups()) {
        if (sys.numTemperatureGroups() != syo.numTemperatureGroups()) {
          cerr << "WARNING: The number of temperature groups has changed. manual check recommended.\n";
          cerr << "         Number of temperature greoups: " << sys.numTemperatureGroups() << " (" << args["topo"] << ")" << endl
                  << "                                        " << syo.numTemperatureGroups() << " (reduced topology)\n";
        }
      }
    }
    
    // and write out the new topology
    OutTopology ot(cout);
    ostringstream os;
    os << "Reduced topology based on " << args["topo"] << endl;
    os << "using atoms ";
    for(unsigned int i=0; i< as_str.size(); i++)
      os << as_str[i] << " ";
    
    ot.setTitle(os.str());
    
    ot.write(syo, it.forceField());
    */
    return 0;
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}



