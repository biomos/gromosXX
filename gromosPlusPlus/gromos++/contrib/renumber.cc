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
 * @file renumber.cc
 * re-assigns force field types in building blocks
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor renumber
 * @section renumber re-assigns force field types in building blocks
 * @author @ref co
 * @date 31-10-08
 *
 * Program renumber can help to convert building blocks from one force field definition to the other. It uses a renumber file, that lists for all bond-, angle-, dihedral-, improper- and atomtypes the original type number and the new type number. Note that it only translate the types for a building block, but does not introduce new types or modify charge distributions.
 *
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@build</td><td>&lt;mtb-file&gt; </td></tr>
 * <tr><td> \@block</td><td>&lt;buildingblock name&gt; </td></tr>
 * <tr><td> \@renum</td><td>&lt;renumber file&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  renumber
    @build   mtb45a4.dat
    @block   ALA
    @renum   ren_45a4_to_53a6.dat
 @endverbatim
 *
 * <hr>
 */
#include <cassert>
#include <cstddef>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/gio/InBuildingBlock.h"
#include "../src/gcore/BuildingBlock.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/BbSolute.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gio/Ginstream.h"
#include "../src/gio/OutBuildingBlock.h"
#include "../src/gromos/Exception.h"

using namespace gcore;
using namespace gio;
using namespace args;
using namespace std;


int main(int argc, char *argv[]){
  Argument_List knowns;
  knowns << "build" << "block" << "renum";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@build <mtb-file>\n";
  usage += "\t@block <buildingblock name>\n";
  usage += "\t@renum <renumber file>\n";

  try{
    Arguments args(argc, argv, knowns, usage);

    // read in the building block file
    InBuildingBlock ibb(args["build"]);
    BuildingBlock mtb(ibb.building());
    
    // define the selected block
    int index=mtb.findBb(args["block"]);
    if(index==0)
      throw gromos::Exception("renumber",
	      "Building block " + args["block"] + " not found" );
    BbSolute bb;
    bool endgroup = false;
    
    if(index>0)
      bb=mtb.bb(index-1);
    else{
      
      bb=mtb.be(-index-1);
      endgroup = true;
    }
    
    // read in the renumber file
    Ginstream gin(args["renum"]);

    // quite ugly to define this here, but who cares
    const int max_number_of_types=100;
    
    map<int, int> atomtypes, bondtypes, angletypes, impropertypes, 
      dihedraltypes;
    std::vector<std::string> buffer;
    std::vector<std::vector<std::string > > content;
    while(!gin.stream().eof()){
      gin.getblock(buffer);
      if(!gin.stream().eof()){ 
	if(buffer[buffer.size()-1].find("END")!=0)
	  throw gromos::Exception("renumber", "Renumber file " + gin.name() +
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
      if      ((*iter)[0]=="ATOMTYPE"     || (*iter)[0]=="ATOMTYPECONV")
	pointer_to_a_map = &atomtypes;
      else if ((*iter)[0]=="BONDTYPE"     || (*iter)[0]=="BONDTYPECONV")
	pointer_to_a_map = &bondtypes;
      else if ((*iter)[0]=="ANGLETYPE"    || (*iter)[0]=="ANGLETYPECONV")
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
	(*pointer_to_a_map)[--a]=--b;
      }
    }
    // let's fill up all types that are not used with themselves

    //atomtypes
    for(int i=1; i< max_number_of_types; i++) {
      if(atomtypes[i]==0)     atomtypes[i]=i;
      if(bondtypes[i]==0)     bondtypes[i]=i;
      if(angletypes[i]==0)    angletypes[i]=i;
      if(impropertypes[i]==0) impropertypes[i]=i;
      if(dihedraltypes[i]==0) dihedraltypes[i]=i;
    }

    // map atom types
    for(int i = 0; i < bb.numAtoms(); ++i) {
      bb.atom(i).setIac(atomtypes[bb.atom(i).iac()]);
    }

    // map rest
    for(BondIterator it(bb); it; ++it) {
      it().setType(bondtypes[it().type()]);
    }
    for(AngleIterator it(bb); it; ++it) {
      it().setType(angletypes[it().type()]);
    }
    for(ImproperIterator it(bb); it; ++it) {
      it().setType(impropertypes[it().type()]);
    }
    for(DihedralIterator it(bb); it; ++it) {
      it().setType(dihedraltypes[it().type()]);
    }

    OutBuildingBlock obb(cout);

    if (endgroup)
      obb.writeSingle(bb, OutBuildingBlock::BBTypeEnd);
    else
      obb.writeSingle(bb, OutBuildingBlock::BBTypeSolute);
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}





