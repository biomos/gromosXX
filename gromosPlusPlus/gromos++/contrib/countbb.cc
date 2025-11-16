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
 * @file countbb.cc
 * Counts building blocks in building block file
 */
/**
 * @page contrib Contrib program documentation
 *
 * @anchor countbb 
 * @section countbb counts building blocks in building block file
 * @author @ref co
 * @date 31-10-08
 *
 * The program countbb counts the number of solute, solvent and end-group
 * building blocks in a molecular building block file. Besides the number
 * countbb prints out a list of these building blocks. It further 
 * checks wether there are identical building blocks. Countbb also provides
 * the highest types for atoms, bonds, angles, impropers and dihedrals 
 * encountered in the building block file as well as those types which were 
 * not found in the building block file. 
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0> 
 * <tr><td> \@build</td><td>&lt; molecular building block file &gt; </td></tr>
 */


#include <cassert>
#include <iostream>
#include <iomanip>
#include <map>
#include <set>
#include <string>

#include "../src/args/Arguments.h"
#include "../src/gio/InBuildingBlock.h"
#include "../src/gcore/BuildingBlock.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/BbSolute.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gromos/Exception.h"


using namespace gcore;
using namespace gio;
using namespace args;
using namespace std;


int main(int argc, char *argv[]){

  Argument_List knowns;
  knowns << "build";
  
  string usage = string("# ") + argv[0];
  usage += "\n\t@build <mtb-file>\n";

  try{
    Arguments args(argc, argv, knowns, usage);

    // read in the building block file
    BuildingBlock mtb;
    Arguments::const_iterator iter=args.lower_bound("build"),
      to=args.upper_bound("build");
    for( ; iter!=to ; ++iter){
      InBuildingBlock ibb(iter->second);
      mtb.addBuildingBlock(ibb.building());
    }

    // keep track of the number of atom, bond, angle etc. types required
    int atomtype=0, bondtype=0, angletype=0, dihedraltype=0, impropertype=0;
    string atomname, bondname, anglename, dihedralname, impropername;
    set<int> atoms, bonds, angles, dihedrals, impropers;
    map<string, int> doubleSolute, doubleSolvent, doubleEnd;
    
    // print out a list of the buildingblocks

    cout << "Number of Solute Building Blocks (MTBUILDBLSOLUTE): "
	 << mtb.numBbSolutes() << endl;
    for(int i=0; i<mtb.numBbSolutes(); i++){
      cout << "  " << mtb.bb(i).resName() << endl;
      // check if it was there before
      for(int j=0; j<i; ++j){
	if(mtb.bb(i).resName()==mtb.bb(j).resName())
	  doubleSolute[mtb.bb(i).resName()]++;
      }
      for(int j=0; j<mtb.bb(i).numAtoms(); j++){
	atoms.insert(mtb.bb(i).atom(j).iac());
	if(mtb.bb(i).atom(j).iac() > atomtype) {
	  atomname=mtb.bb(i).resName();
	  atomtype=mtb.bb(i).atom(j).iac();
	}
      }
      
      BondIterator bi(mtb.bb(i));
      for(;bi;++bi){
	bonds.insert(bi().type());
	if(bi().type()>bondtype){
	  bondname= mtb.bb(i).resName();
	  bondtype=bi().type();
	}
      }
      
      AngleIterator ai(mtb.bb(i));
      for(; ai; ++ai){
	angles.insert(ai().type());
	if(ai().type()>angletype){
	  anglename=mtb.bb(i).resName();
	  angletype=ai().type();
	}
      }
      
      ImproperIterator ii(mtb.bb(i));
      for(; ii; ++ii){
	impropers.insert(ii().type());
	if(ii().type()>impropertype){
	  impropername=mtb.bb(i).resName();
	  impropertype=ii().type();
	}
      }
      
      DihedralIterator di(mtb.bb(i));
      for(; di; ++di){
	dihedrals.insert(di().type());
	if(di().type()>dihedraltype){
	  dihedralname=mtb.bb(i).resName();
	  dihedraltype=di().type();
	}
      }
    }
    cout << endl;
    if(doubleSolute.size()){
      cout << "WARNING the following Solute Building Block";
      if(doubleSolute.size()>1) cout << "s";
      cout << " appear";
      if(doubleSolute.size()==1) cout << "s";
      cout << " more than once:\n";
      map<string,int>::const_iterator iter=doubleSolute.begin(),
	to=doubleSolute.end();
      for(; iter!=to; ++iter){
	cout << "  " << iter->first << " (" << iter->second+1 << ")\n";
      }
      cout << endl;
    }
    
    cout << "Number of Solvent Building Blocks (MTBUILDBLSOLVENT): "
	 << mtb.numBbSolvents() << endl;
    for(int i=0; i<mtb.numBbSolvents(); i++){
      cout << "  " << mtb.bs(i).solvName() << endl;
      // check if it was there before
      for(int j=0; j<i; ++j){
	if(mtb.bs(i).solvName()==mtb.bs(j).solvName())
	  doubleSolvent[mtb.bs(i).solvName()]++;
      }
      for(int j=0; j<mtb.bs(i).numAtoms(); j++){
	atoms.insert(mtb.bs(i).atom(j).iac());
	if(mtb.bs(i).atom(j).iac() > atomtype){
	  atomname=mtb.bs(i).solvName();
	  atomtype=mtb.bs(i).atom(j).iac();
	}
      }
    }
    cout << endl;
    if(doubleSolvent.size()){
      cout << "WARNING the following Solvent Building Block";
      if(doubleSolvent.size()>1) cout << "s";
      cout << " appear";
      if(doubleSolvent.size()==1) cout << "s";
      cout << " more than once:\n";
      map<string,int>::const_iterator iter=doubleSolvent.begin(),
	to=doubleSolvent.end();
      for(; iter!=to; ++iter){
	cout << "  " << iter->first << " (" << iter->second+1 << ")\n";
      }
      cout << endl;
    }    
    
    cout << "Number of End-group Building Blocks (MTBUILDBLEND): "
	 << mtb.numBbEnds() << endl;
    for(int i=0; i<mtb.numBbEnds(); i++){
	
      cout << "  " << mtb.be(i).resName() << endl;
      // check if it was there before
      for(int j=0; j<i; ++j){
	if(mtb.be(i).resName()==mtb.be(j).resName())
	  doubleEnd[mtb.be(i).resName()]++;
      }
      for(int j=0; j<mtb.be(i).numAtoms(); j++){
	
	atoms.insert(mtb.be(i).atom(j).iac());
	if(mtb.be(i).atom(j).iac() > atomtype) {
	  atomname=mtb.be(i).resName();
	  atomtype=mtb.be(i).atom(j).iac();
	}
      }
      
      BondIterator bi(mtb.be(i));
      for(;bi;++bi){
	bonds.insert(bi().type());
	if(bi().type()>bondtype){
	  bondname=mtb.be(i).resName();
	  bondtype=bi().type();
	}
      }
      
      AngleIterator ai(mtb.be(i));
      for(; ai; ++ai){
	angles.insert(ai().type());
	if(ai().type()>angletype){
	  anglename=mtb.be(i).resName();
	  angletype=ai().type();
	}
      }
      
      ImproperIterator ii(mtb.be(i));
      for(; ii; ++ii){
	impropers.insert(ii().type());
	if(ii().type()>impropertype){
	  impropername=mtb.be(i).resName();
	  impropertype=ii().type();
	}
      }
      
      DihedralIterator di(mtb.be(i));
      for(; di; ++di){
	dihedrals.insert(di().type());
	if(di().type()>dihedraltype){
	  dihedralname=mtb.be(i).resName();
	  dihedraltype=di().type();
	}
      }
    }
    cout << endl;

    if(doubleEnd.size()){
      cout << "WARNING the following End-group Building Block";
      if(doubleEnd.size()>1) cout << "s";
      cout << " appear";
      if(doubleEnd.size()==1) cout << "s";
      cout << " more than once:\n";
      map<string,int>::const_iterator iter=doubleEnd.begin(),
	to=doubleEnd.end();
      for(; iter!=to; ++iter){
	cout << "  " << iter->first << " (" << iter->second+1 << ")\n";
      }
      cout << endl;
    }    

    cout << "Highest types encountered for:   (in building block)" << endl;
    cout << "  atoms     : " << setw(5) << atomtype+1 
	 << "  (" << atomname << ")" << endl;
    cout << "  bonds     : " << setw(5) << bondtype +1 
	 << "  (" << bondname << ")" << endl;
    cout << "  angles    : " << setw(5) << angletype +1 
	 << "  (" << anglename << ")" << endl;
    cout << "  impropers : " << setw(5) << impropertype +1 
	 << "  (" << impropername << ")"<< endl;
    cout << "  dihedrals : " << setw(5) << dihedraltype +1 
	 << "  (" << dihedralname << ")" << endl;
    cout << endl;
    cout << "Types that were not seen below these values: " << endl;
    cout << "  atoms     : ";
    for(int i=0; i<atomtype; i++)
      if(atoms.count(i)==0) cout << setw(5) << i+1 << " ";
    cout << endl;
    
    cout << "  bonds     : ";
    for(int i=0; i<bondtype; i++)
      if(bonds.count(i)==0) cout << setw(5) <<i+1 << " ";
    cout << endl;
    
    cout << "  angles    : ";
    for(int i=0; i<angletype; i++)
      if(angles.count(i)==0) cout << setw(5) <<i+1 << " ";
    cout << endl;
    
    cout << "  impropers : ";
    for(int i=0; i<impropertype; i++)
      if(impropers.count(i)==0) cout << setw(5) <<i+1 << " ";
    cout << endl;
    
    cout << "  dihedrals : ";
    for(int i=0; i<dihedraltype; i++)
      if(dihedrals.count(i)==0) cout << setw(5) <<i+1 << " ";
    cout << endl;
    
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}





