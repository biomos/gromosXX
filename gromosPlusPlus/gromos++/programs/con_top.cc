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
 * @file con_top.cc
 * Convert topology to different force field version
 */

/**
 * @page programs Program Documentation
 *
 * @anchor con_top
 * @section con_top Convert topology to different force field version
 * @author @ref co
 * @date 6-6-07
 *
 * A molecular topology file in which a system is described by a specific 
 * version of a force field parameter set can be converted into a molecular 
 * topology file with interaction parameters from a different force field 
 * version, using the program con_top.
 *
 * An interaction function parameter file has to be specified that correspond 
 * to the force field version into which the molecular topology should be 
 * converted. con_top checks whether the topology is not referring to atom, 
 * bond, bond angle or (improper) dihedral types that are not defined in the 
 * parameter file. If type numbers of atoms, bonds, bond angles, etc. change 
 * with the force field parameter set, a renumbering file (see volume 4 for 
 * its formats) can be given to specify these changes.
 * 
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file to be converted&gt; </td></tr>
 * <tr><td> \@param</td><td>&lt;interaction function parameter file&gt; </td></tr>
 * <tr><td> \@renum</td><td>&lt;renumbering file&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  con_top
    @topo   ex.top
    @param  ifp53a6.dat
    @renum  ren45a4_to_53a5.dat
 @endverbatim
 *
 * <hr>
 */

#include <cstddef>
#include <iostream>
#include <cassert>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InParameter.h"
#include "../src/gio/OutTopology.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/CrossDihedral.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gio/Ginstream.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;

void renumber_types(System &sys, GromosForceField &gff, string renum);
void check_types(System &sys, GromosForceField &gff);

int main(int argc, char *argv[]){

  Argument_List knowns;
  knowns << "topo" << "param" << "renum";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo  <molecular topology file to be converted>\n";
  usage += "\t@param <interaction function parameter file>\n";
  usage += "\t@renum <renumbering file>\n";
  
  try{
    Arguments args(argc, argv, knowns, usage);

    InTopology it(args["topo"]);

    System sys(it.system());
    
    OutTopology ot(cout);
    string addtitle;
    addtitle+="CON_TOP parameters: "+args["param"];
    if (args.count("renum")>0) 
      addtitle += "\nrenumbering from: " + args["renum"];
    
    ot.setTitle(it.title()+addtitle);

    InParameter ip(args["param"]);
    GromosForceField gff(ip.forceField());
    gff.setFpepsi(it.forceField().fpepsi());
    gff.setHbar(it.forceField().hbar());
    gff.setSpdl(it.forceField().spdl());
    gff.setBoltz(it.forceField().boltz());

    // maybe some types are to be renumbered?
    if(args.count("renum")>0)
      renumber_types(sys, gff, args["renum"]);

    // check if the topology is now not referring to non-existing types
    check_types(sys, gff);
    
    ot.write(sys,gff);
    
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}

void check_types(System &sys, GromosForceField &gff)
{
  int num_atom_types=gff.numAtomTypeNames();
  int num_bond_types=gff.numBondTypes();
  int num_angle_types=gff.numAngleTypes();
  int num_improper_types=gff.numImproperTypes();
  int num_dihedral_types=gff.numDihedralTypes();
  
  for(int m=0; m< sys.numMolecules();m++){
    for(int a=0; a<sys.mol(m).numAtoms(); a++){
      if(sys.mol(m).topology().atom(a).iac() >= num_atom_types){
	ostringstream os;
	os << "Atom " << m+1 << ":" << a+1 << " has a higher integer "
	   << "atom code (" << sys.mol(m).topology().atom(a).iac()+1
	   << ") than defined in the parameter file ("
	   << num_atom_types << ")";
	
	throw gromos::Exception("con_top", os.str());
      }
    }
    BondIterator bi(sys.mol(m).topology());
    for(;bi;++bi){
      if(bi().type() >= num_bond_types){
	ostringstream os;
	os << "Bond " << m+1 << ":" << bi()[0] << " - " << m+1 << ":" 
	   << bi()[1] << " has a higher bondtype (" << bi().type()+1 
	   << ") than defined in the parameter file ("
	   << num_bond_types << ")";
	throw gromos::Exception("con_top", os.str());
      }
    }
    AngleIterator ai(sys.mol(m).topology());
    for(;ai;++ai){
      if(ai().type() >= num_angle_types){
	ostringstream os;
	os << "Angle " << m+1 << ":" << ai()[0]+1 << " - " << m+1 << ":" 
	   << ai()[1]+1 << " - " << m+1 << ":" << ai()[2]+1
	   << " has a higher angletype (" << ai().type()+1 
	   << ") than defined in the parameter file ("
	   << num_angle_types << ")";
	throw gromos::Exception("con_top", os.str());
      }
    }
    ImproperIterator ii(sys.mol(m).topology());
    for(;ii;++ii){
      if(ii().type() >= num_improper_types){
	ostringstream os;
	os << "Improper " << m+1 << ":" << ii()[0]+1 << " - " << m+1 << ":" 
	   << ii()[1]+1 << " - " << m+1 << ":" << ii()[2]+1 << " - " 
	   << m+1 << ":" << ii()[3]+1 
	   << " has a higher impropertype (" << ii().type()+1 
	   << ") than defined in the parameter file ("
	   << num_improper_types << ")";
	throw gromos::Exception("con_top", os.str());
      }
    }    
    DihedralIterator di(sys.mol(m).topology());
    for(;di;++di){
      if(di().type() >= num_dihedral_types){
	ostringstream os;
	os << "Dihedral " << m+1 << ":" << di()[0]+1 << " - " << m+1 << ":"
	   << di()[1]+1 << " - " << m+1 << ":" << di()[2]+1 << " - " 
	   << m+1 << ":" << di()[3]+1 
	   << " has a higher dihedraltype (" << di().type()+1 
	   << ") than defined in the parameter file ("
	   << num_improper_types << ")";
	throw gromos::Exception("con_top", os.str());
      }
    }

    CrossDihedralIterator cdi(sys.mol(m).topology());
    for(;cdi;++cdi){
      if(cdi().type() >= num_dihedral_types){
	ostringstream os;
	os << "CrossDihedral " << m+1 << ":" << cdi()[0]+1 << " - " << m+1 << ":"
	   << cdi()[1]+1 << " - " << m+1 << ":" << cdi()[2]+1 << " - "
	   << m+1 << ":" << cdi()[3]+1
	   << " has a higher dihedraltype (" << cdi().type()+1
	   << ") than defined in the parameter file ("
	   << num_improper_types << ")";
	throw gromos::Exception("con_top", os.str());
      }
    }
  }
}

void renumber_types(System &sys, GromosForceField &gff, string renum)
{
  // read in the renumber file
  Ginstream gin(renum);

  int num_atom_types=gff.numAtomTypeNames();
  int num_bond_types=gff.numBondTypes();
  int num_angle_types=gff.numAngleTypes();
  int num_improper_types=gff.numImproperTypes();
  int num_dihedral_types=gff.numDihedralTypes();
    
  map<int, int> atomtypes, bondtypes, angletypes, impropertypes, 
    dihedraltypes;
  std::vector<std::string> buffer;
  std::vector<std::vector<std::string > > content;
  while(!gin.stream().eof()){
    gin.getblock(buffer);
    if(!gin.stream().eof()){
      if(buffer[buffer.size()-1].find("END")!=0)
	throw gromos::Exception("con_top", "Renumber file " + gin.name() +
				" is corrupted. No END in "+buffer[0]+
				" block. Got\n"
				+ buffer[buffer.size()-1]);
      content.push_back(buffer);
    }    
  }
  // now loop over the content
  std::vector<std::vector<std::string > >::const_iterator 
    iter=content.begin();
  for( ; iter!=content.end(); ++iter){
    map<int, int> *pointer_to_a_map=NULL;
    if      ((*iter)[0]=="ATOMTYPE"    || (*iter)[0]=="ATOMTYPECONV") 
      pointer_to_a_map = &atomtypes;
    else if ((*iter)[0]=="BONDTYPE"    || (*iter)[0]=="BONDTYPECONV")
      pointer_to_a_map = &bondtypes;
    else if ((*iter)[0]=="ANGLETYPE"   || (*iter)[0]=="ANGLETYPECONV")
      pointer_to_a_map = &angletypes;
    else if ((*iter)[0]=="IMPROPERTYPE" || (*iter)[0]=="IMPROPERTYPECONV")
      pointer_to_a_map = &impropertypes;
    else if ((*iter)[0]=="DIHEDRALTYPE" || (*iter)[0]=="DIHEDRALTYPECONV")
      pointer_to_a_map = &dihedraltypes;
    else
      throw gromos::Exception("renumber", 
			      "Don't know how to handle "+(*iter)[0] + "-block");
    
    int a, b;
    
    // now loop over the contents of the block
    for(unsigned int i=1; i< (*iter).size()-1; i++){
      std::istringstream linestream((*iter)[i]);
      linestream >> a >> b;
      (*pointer_to_a_map)[a]=b;
    }
  }
  // let's fill up all types that are not used with themselves
  
  //atomtypes
  for(int i=1; i< num_atom_types; i++)
    if(atomtypes[i]==0) atomtypes[i]=i;
  for(int i=1; i< num_bond_types; i++) 
    if(bondtypes[i]==0) bondtypes[i]=i;
  for(int i=1; i< num_angle_types; i++) 
    if(angletypes[i]==0) angletypes[i]=i;
  for(int i=1; i< num_improper_types; i++) 
    if(impropertypes[i]==0) impropertypes[i]=i;
  for(int i=1; i< num_dihedral_types; i++) 
    if(dihedraltypes[i]==0) dihedraltypes[i]=i;
  
  // Now loop over all the bonds, angles and atoms in the topology to
  // replace the types

  // molecules
  for(int m=0; m<sys.numMolecules(); m++){

    MoleculeTopology mt;
    
    // atoms
    for(int a=0; a < sys.mol(m).numAtoms(); a++){
      
      sys.mol(m).topology().atom(a).setIac( 
	atomtypes[sys.mol(m).topology().atom(a).iac()+1] - 1);
      mt.addAtom(sys.mol(m).topology().atom(a));
      mt.setResNum(a,sys.mol(m).topology().resNum(a));
    }
    for(int r=0; r < sys.mol(m).topology().numRes(); r++)
      mt.setResName(r, sys.mol(m).topology().resName(r));
    
    // bonds
    BondIterator bi(sys.mol(m).topology());
    for(;bi;++bi){
      Bond b=bi();
      b.setType(bondtypes[bi().type()+1] - 1);
      mt.addBond(b);
    }
    
    // Angles
    AngleIterator ai(sys.mol(m).topology());
    for(;ai; ++ai){
      Angle a=ai();
      a.setType(angletypes[ai().type()+1] - 1);
      mt.addAngle(a);
    }
    
    // Impropers
    ImproperIterator ii(sys.mol(m).topology());
    for(;ii; ++ii){
      Improper i=ii();
      i.setType(impropertypes[ii().type()+1] - 1);
      mt.addImproper(i);
    }
    
    // Dihedrals
    DihedralIterator di(sys.mol(m).topology());
    for(;di; ++di){
      Dihedral d=di();
      d.setType(dihedraltypes[di().type()+1] - 1);
      mt.addDihedral(d);
    }

    // CrossDihedrals
    CrossDihedralIterator cdi(sys.mol(m).topology());
    for(;cdi; ++cdi){
      CrossDihedral d=cdi();
      d.setType(dihedraltypes[cdi().type()+1] - 1);
      mt.addCrossDihedral(d);
    }
    sys.mol(m).topology() = mt;

    // after possibly renumbering atom types, or new masses
    // we have to assign hydrogen atoms new.
    sys.mol(m).topology().clearH();
    // determine H atoms based on the masses
    sys.mol(m).topology().setHmass(1.008);
    // or on the iac-code
    // sys.mol(m).topology().setHmass(17);
  }
  // And don't forget the solvent!! Okay, I did forget the solven.
  for(int s=0; s<sys.numSolvents(); s++){
    for(int a=0; a<sys.sol(s).topology().numAtoms(); a++){
      sys.sol(s).topology().atom(a).setIac(
	   atomtypes[sys.sol(s).topology().atom(a).iac()+1] - 1);
    }
  }
  
}








