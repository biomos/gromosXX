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
 * @file make_top.cc
 * Creates molecular topologies from building blocks
 */

/**
 * @page programs Program Documentation
 *
 * @anchor make_top
 * @section make_top Creates molecule topologies from building blocks
 * @author @ref co
 * @date 5.6.07
 *
 * Program make_top builds a molecular topology from a building block sequence.
 * make_top reads in a molecular topology building-block file (e.g. 
 * mtb53a6.dat) and an interaction function parameter file (e.g. ifp53a6.dat), 
 * and gathers the specified building blocks to create a topology. Cysteine 
 * residues involved in disulfide bridges as well as heme and coordinating 
 * residues involved in covalent bonds to the iron atom have to be explicitly 
 * specified. Topologies for cyclic sequences of building blocks can be 
 * generated using the keyword 'cyclic'. The resulting topology file is written
 * out to the standard output.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@build</td><td>&lt;molecular topology building block file&gt; </td></tr>
 * <tr><td> \@param</td><td>&lt;interaction function parameter file&gt; </td></tr>
 * <tr><td> \@seq</td><td>&lt;sequence of building blocks in the solute&gt; </td></tr>
 * <tr><td> \@solv</td><td>&lt;building block for the solvent&gt; </td></tr>
 * <tr><td> [\@cys</td><td>&lt;cys1&gt;-&lt;cys2&gt; .. &lt;cys1&gt;-&lt;cys2&gt;] </td></tr>
 * <tr><td> [\@heme</td><td>&lt;residue sequence number&gt; &lt;heme sequence number&gt;] </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  make_top
    @build    mtb53a6.dat
    @param    ifp53a6.dat
    @seq      NH3+ ALA CYS1 GLU HIS1 CYS2 GLY COO- HEME NA+
    @solv     H2O
    @cys      2-5
    @heme     4 7
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <map>
#include <set>
#include <cstdio>

#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InParameter.h"
#include "../src/gio/InBuildingBlock.h"
#include "../src/gio/OutTopology.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/BuildingBlock.h"
#include "../src/gcore/BbSolute.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/LinearTopology.h"
#include "../src/gromos/Exception.h"
#include "../src/utils/make_top.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;

int main(int argc, char *argv[]){

  Argument_List knowns;
  knowns << "build" << "param" << "seq" << "solv" << "cys" << "heme";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@build <molecular topology building block file>\n";
  usage += "\t@param <interaction function parameter file>\n";
  usage += "\t@seq   <sequence of building blocks in the solute>\n";
  usage += "\t@solv  <building block for the solvent>\n";
  usage += "\t[@cys  <cys1>-<cys2> .. <cys1>-<cys2>]\n";
  usage += "\t       (sequence numbers for disulfide bridges)\n";
  usage += "\t[@heme <residue sequence number> <heme sequence number>]\n";
  usage += "\t       (to covalently bind a heme group)\n";
  
  try{
    Arguments args(argc, argv, knowns, usage);

    // read in the force field parameter file
    InParameter ip(args["param"]);
    GromosForceField gff(ip.forceField());

    // read in the building block file
    BuildingBlock mtb;
    Arguments::const_iterator iter=args.lower_bound("build"),
      to=args.upper_bound("build");
    for( ; iter!=to ; ++iter){
      InBuildingBlock ibb(iter->second);
      mtb.addBuildingBlock(ibb.building());
      
      // Check force field consistency
      std::set<std::string> FF_set = std::move(ibb.building().ForceField());
      if (FF_set.find(gff.ForceField()) == FF_set.end()){
        std::string FF_codes;
        std::set<std::string>::iterator itr;
        for (itr = FF_set.begin(); itr != FF_set.end(); itr++){
          FF_codes += *itr + " ";
        }
        throw gromos::Exception("make_top", 
			      "Parameter file and building block file(s) have "
			      "different FORCEFIELD codes\nParameter file: "
			      +gff.ForceField()
			      + "\nBuilding block file: " + FF_codes);
      }
    }
    
    // parse the input for disulfide bridges    
    vector<int> cys1, cys2, csa1, csa2;
 
    for(Arguments::const_iterator iter=args.lower_bound("cys"),
	  to=args.upper_bound("cys"); iter!=to; ++iter){
      string s=iter->second;
      
      int a;
      std::string::size_type iterator;
      iterator = s.find('-');
      if (iterator == std::string::npos)
        throw gromos::Exception("make_top", "Bad cysteine specification\n");
      if (sscanf((s.substr(0,iterator)).c_str(),"%d", &a) != 1)
        throw gromos::Exception("make_top", 
	       "Bad first cysteine specification: "+s+"\n");
      cys1.push_back(--a);
      if (sscanf((s.substr(iterator+1,s.length())).c_str(),"%d", &a) != 1)
        throw gromos::Exception("make_top", 
		"Bad second cysteine specification: "+s+"\n");
      cys2.push_back(--a);
    }
 
    // parse the input for heme groups
    vector<int> his1, heme, hsa1, hsn2, hma;
    vector<string> atomToHeme;
    
    for(Arguments::const_iterator iter=args.lower_bound("heme"),
	  to = args.upper_bound("heme"); iter != to; ++iter){
      int h1=atoi(iter->second.c_str());
      ++iter;
      if(iter==to)
	throw gromos::Exception("make_top", 
				"Bad heme-linking specification: give pairs");
      int h2=atoi(iter->second.c_str());
      his1.push_back(--h1);
      heme.push_back(--h2);
    }
    
    //some variables and lists to store data temporarily
    int index=0;
    //status=0: normal linking
    //status=1: current is a beginning
    //status=2: current is first after beginning
    //status=3: current is an end
    int status=0;
    int repforward=0;
    // firstAtom is the first atom of the current molecule as determined
    // from a starting end-group
    int firstAtom=0;
    
    gcore::LinearTopology lt;
    int resnum=0;
    int cyclic=0;
  
    
    // loop over the sequence
    for(Arguments::const_iterator iter=args.lower_bound("seq"),
	  to=args.upper_bound("seq"); iter!=to; ++iter){
      
      if(iter->second == "cyclic"){
	if(lt.atoms().size())
	  throw(gromos::Exception("make_top", 
				  "make_top can only cyclize one complete "
				  "molecule. The keyword cyclic should be the "
				  "first in the sequence"));
        prepareCyclization(lt);
        iter++;
	status = 1;
        repforward = 0;
	cyclic=1;
      }
      int countBB = 0;
      index = mtb.findBb(iter->second, countBB);
      if(index==0) throw gromos::Exception("make_top", 
					   "Cannot find building block for "
					   +iter->second+
					   " in building block file(s)");
      if(countBB!=1) 
	cerr << "WARNING: Found more than one version of building block for "
	     << iter->second << ".\n"
	     << "Using the first that was encountered.\n\n";
      
      //determine index and status
      if(index<0) {
	index=-1-index;
        if(mtb.be(index).rep()<0) status=3;
	else status = 1;
      }
      else {
	index=index-1;
	if (status==1) status =2;
	else status = 0;
      }

      // depending on the status add the correct building block to the
      // linearTopology
      switch(status){
      case 0:
        addSolute(lt, mtb.bb(index), resnum, iter->second, 0, firstAtom);
        resnum++;
	break;
      case 1:
        repforward = addBegin(lt, mtb.be(index), resnum);
	firstAtom = lt.atoms().size() - mtb.be(index).numAtoms();
        addCovEnd(lt, mtb.be(index), firstAtom);
	break;
      case 2:
        addSolute(lt, mtb.bb(index), resnum, iter->second,
		  repforward, firstAtom);
	// a call to remove atoms, because we might have some negative iac's
	// from the beginning buildingblock.
        lt.removeAtoms();
	resnum++;
	break;
      case 3:
	// this residue actually belongs to the previous one
        //check if there is a previous one
        if ( resnum ==0 ){
           throw gromos::Exception("make_top", 
			      "End group " 
                              +iter->second+ 
                              " cannot be at start of sequency");
        }
        resnum--;
        addEnd(lt, mtb.be(index), resnum);
        addCovEnd(lt,mtb.be(index),lt.atoms().size()-mtb.be(index).numAtoms());
	resnum++;
	break;
      }
    }
          
    // this would be the place to handle any cysteine bridges
    for(unsigned int j=0; j<cys1.size(); j++){
      bool found=false;
      for(unsigned int k=0; k<lt.resMap().size();k++){
	if(lt.resMap()[k]==cys1[j]&&lt.atoms()[k].name()=="CA") 
	  {csa1.push_back(k); found=true;
          if(lt.resNames()[lt.resMap()[k]]!="CYS1") {
          cerr << "WARNING: residue number " << int(cys1[j]+1) << " is not named CYS1. \n"
               << "Residue name: " << lt.resNames()[lt.resMap()[k]] << "\n";
           }
       break;}}
     
      if(found==false){
      ostringstream os;
        os << "Could not find a residue numbered: "
           << int(cys1[j]+1) << "\n"
           << "as specified by @cys. \n";
        throw gromos::Exception("make_top",os.str());
      }
    }

    for(unsigned int j=0; j<cys2.size(); j++){
      bool found=false;
//      cout << "this is cys2 number" << cys2[j] << endl;
      for(unsigned int k=0; k<lt.resMap().size();k++){
	if(lt.resMap()[k]==cys2[j]&&lt.atoms()[k].name()=="CA") 
	  {csa2.push_back(k); found=true;
          if(lt.resNames()[lt.resMap()[k]]!="CYS2") {
        cerr << "WARNING: residue number " << int(cys2[j]+1) << " is not named CYS2. \n"
             << "Residue name: " << lt.resNames()[lt.resMap()[k]] << "\n";
         }
       break;}}

      if(found==false){
      ostringstream os;
        os << "Could not find a residue numbered: "
           << int(cys2[j]+1) << "\n"
           << "as specified by @cys. \n";
      throw gromos::Exception("make_top",os.str());
     }
    }
    

    for(unsigned int j=0; j<csa1.size(); j++)
      setCysteines(lt, csa1[j], csa2[j]);
    
    // and maybe an irritating heme group?
    for(unsigned int j=0; j<his1.size(); j++)
      for(unsigned int k=0; k<lt.resMap().size(); k++)
	if(lt.resMap()[k]==his1[j] && lt.atoms()[k].name()=="CA")
	  { hsa1.push_back(k); break;}
    for(unsigned int j=0; j<his1.size(); j++)
      for(unsigned int k=0; k<lt.resMap().size(); k++)
	if(lt.resMap()[k]==his1[j] && 
	   (lt.atoms()[k].name()=="NE2" || lt.atoms()[k].name()=="SG"))
	  { hsn2.push_back(k); break;}
    for(unsigned int j=0; j<heme.size(); j++)
      for(unsigned int k=0; k<lt.resMap().size(); k++)
	if(lt.resMap()[k]==heme[j])
	  { hma.push_back(k); break;}
    if(hsa1.size() != his1.size() || hsn2.size() != his1.size())
      throw gromos::Exception("make_top", 
			      "Residues to connect to heme requires an atom CA "
			      "and an atom NE2 / SG. One of these is not "
			      "found.");
  
    if(hma.size() != his1.size())
      throw gromos::Exception("make_top", 
			      "For covalent interaction to Heme, an atom "
			      "called Fe is required");
    
    for(unsigned int j=0; j<hsa1.size(); j++)
      setHeme(lt, hsa1[j], hsn2[j], hma[j]);
	    
    // possibly cyclize
    if(cyclic) cyclize(lt);
    
    // get the 1-4 interactions from the bonds
    lt.get14s();
  
    // transform masses from integer to double
    for(unsigned int i=0; i< lt.atoms().size(); i++){
      double m=gff.findMass(int(lt.atoms()[i].mass()));
      if(m!=0.0)
	lt.atoms()[i].setMass(m);
      else{
	ostringstream os;
	os << "Could not find masstype " 
	   << int(lt.atoms()[i].mass()+1) 
	   << " in parameter file (atom " << i+1 << "; "
	   << lt.atoms()[i].name() << ").";
	throw gromos::Exception("make_top",os.str());
      }
    }
    
    // do some checks before preparing the linear topology to be written out
    //
    int numAtoms = lt.atoms().size();
    // do the bonds make sense?
    //
    // no checks neede for the bonds here, they are already checked in linearTopology
    //
    // are the exclusions within the solute?
    for (unsigned int a = 0; a < lt.atoms().size(); ++a) {
      for (int e = 0; e < lt.atoms()[a].exclusion().size(); ) {
        if(lt.atoms()[a].exclusion().atom(e) >= (int)lt.atoms().size()) {
          cerr << "WARNING: exclusion SKIPPED since it is not within the solute:\n";
          cerr << "         atoms " << a + 1 << " had atom " << 
                  lt.atoms()[a].exclusion().atom(e) + 1 << " in its exclusion list\n";
          lt.atoms()[a].exclusion().erase(lt.atoms()[a].exclusion().atom(e));
        } else {
          ++e;
        }
      }
    }
    //
    // do the bond angle make sense?
    for(set<Angle>::const_iterator it = lt.angles().begin();
            it != lt.angles().end();) {
      if(((*it)[0] >= numAtoms) || ((*it)[1] >= numAtoms) || ((*it)[2] >= numAtoms)
              || (*it)[0] < 0 || (*it)[1] < 0 || (*it)[2] < 0) {
        cerr << "WARNING: bond angle SKIPPED since it is not within the solute:\n";
        cerr << "         " << (*it)[0] + 1 << "-" << (*it)[1] + 1 << "-" << (*it)[2] + 1 << endl;
        it = lt.angles().erase(it);
      } else {
        ++it;
      }
    }
    // do the improper dihedral make sense?
    for (set<Improper>::const_iterator it = lt.impropers().begin();
            it != lt.impropers().end();) {
      if(((*it)[0] >= numAtoms) || ((*it)[1] >= numAtoms) || ((*it)[2] >= numAtoms) || ((*it)[3] >= numAtoms)
              || (*it)[0] < 0 || (*it)[1] < 0 || (*it)[2] < 0 || (*it)[3] < 0){
        cerr << "WARNING: improper dihedral SKIPPED since it is not within the solute:\n";
        cerr << "         " << (*it)[0] + 1 << "-" << (*it)[1] + 1 << "-" << (*it)[2] + 1 << "-" << (*it)[3] + 1 << endl;
        it = lt.impropers().erase(it);
      } else {
        ++it;
      }
    }
    // do the dihedral make sense?
    for (set<Dihedral>::const_iterator it = lt.dihedrals().begin();
            it != lt.dihedrals().end(); ) {
      if(((*it)[0] >= numAtoms) || ((*it)[1] >= numAtoms) || ((*it)[2] >= numAtoms) || ((*it)[3] >= numAtoms)
              || (*it)[0] < 0 || (*it)[1] < 0 || (*it)[2] < 0 || (*it)[3] < 0){
        cerr << "WARNING: dihedral SKIPPED since it is not within the solute:\n";
        cerr << "         " << (*it)[0] + 1 << "-" << (*it)[1] + 1 << "-" << (*it)[2] + 1 << "-" << (*it)[3] + 1 << endl;
        it = lt.dihedrals().erase(it);
      } else {
        ++it;
      }
    }
    // in case of carbo force-field files, we perform the deletion of dihedral angles with negative types
    set<Dihedral>::const_iterator deleter = lt.dihedrals().begin();
    while(deleter != lt.dihedrals().end()) {
      // in case there is a deleter, go and look for the dihedral to be deleted
      if (deleter->type() < 0) {
        int a1 = (*deleter)[0];
        int a2 = (*deleter)[1];
        int a3 = (*deleter)[2];
        int a4 = (*deleter)[3];
        int type = deleter->type() + 2; // since GROMOS counting starts at 0, usually we ust need to subtract 1, but in the
                                        // case of negative dihedral types we should add on, so add 2
        set<Dihedral>::iterator deletion = lt.dihedrals().begin();
        while(deletion != lt.dihedrals().end()) {
          if ((*deletion)[0] == a1 && (*deletion)[1] == a2 && (*deletion)[2] == a3 && (*deletion)[3] == a4 && deletion->type() == -type) {
            lt.dihedrals().erase(deletion++);
          } else {
            deletion++;
          }
        }
        lt.dihedrals().erase(deleter++);
      } else {
        deleter++;
      }
    }
    
    // do the LJ exceptions make sense?
    for(set<LJException>::const_iterator it = lt.ljexceptions().begin();
            it != lt.ljexceptions().end(); it++) {
      if(((*it)[0] >= numAtoms) || ((*it)[1] >= numAtoms) || (*it)[0] < 0 || (*it)[1] < 0) {
        cerr << "WARNING: LJ exception SKIPPED since it is not within the solute:\n";
        lt.ljexceptions().erase(it);
        cerr << "         " << (*it)[0] + 1 << "-" << (*it)[1] + 1 << endl;
      } 
      // is the condition for the LJ exception fulfilled?
      int numCond = it->numcond();
      if (numCond > 0) {
        bool a1 = false; // the condition or atom 1 is not fulfilled yet...
        bool a2 = false; // the condition or atom 2 is not fulfilled yet...
        int iac1 = lt.atoms()[(*it)[0]].iac();
        int iac2 = lt.atoms()[(*it)[1]].iac();
        for(set<int>::iterator sit = it->cond().begin(); sit != it->cond().end(); sit++) {
          if(*sit == iac1) {
            a1 = true;
          }
          if(*sit == iac2) {
            a2 = true;
          }
          //cout << " cond: " << *sit;
        }
        //cout << " flag: " << it->indicate();
        //cout << " a1: " << (*it)[0] << " " << iac1 << " " << a1;
        //cout << " a2: " << (*it)[1] << " " << iac2 << " " << a2 << " "; 
        // remove if conditions not fulfilled
        if ((!a1 || !a2) && it->indicate() == 0) {
          set<LJException>::const_iterator de=it;
          it--;
          lt.ljexceptions().erase(de);
          cerr << "Removed\n" ;
          cerr << "         " << (*de)[0] + 1 << "-" << (*de)[1] + 1 << endl;
          //it--;
        } else if (!a1 && it->indicate() == 1) {
          set<LJException>::const_iterator de=it;
          it--;
          lt.ljexceptions().erase(de);
          cerr << "Removed\n";
          cerr << "         " << (*de)[0] + 1 << "-" << (*de)[1] + 1 << endl;
          //it--;
        } else if (!a2 && it->indicate() == 2) {
          set<LJException>::const_iterator de=it;
          it--;
          lt.ljexceptions().erase(de);
          cerr << "Removed\n";
          cerr << "         " << (*de)[0] + 1 << "-" << (*de)[1] + 1 << endl;
          //it--;
        }
      }
      // remove the LJ Exceptions from the 14-neighbour list if necessary
      int at1 = (*it)[0];
      int at2 = (*it)[1];
      int n14 = lt.atoms()[at1].exclusion14().size();
      for(int i = 0; i < n14; i++) {
        if(lt.atoms()[at1].exclusion14().atom(i) == at2) {
          lt.atoms()[at1].exclusion14().erase(at2);
        }
      }
    }

    // parse everything into a system    
    System sys;
    lt.parse(sys);
    
    // add the solvent topology
    int countBS = 0;
    index = mtb.findBs(args["solv"], countBS);
    if (index == 0) {
      throw gromos::Exception("make_top",
              "Cannot find building block for "
              + args["solv"] + " in " + args["build"]);
    }
    if (countBS != 1) {
      cerr << "WARNING: Found more than one version of building block for "
              << args["solv"] << ".\n"
              << "Using the first that was encountered.\n\n";
    }
    
    SolventTopology st;

    // adapt the masses
    for(int i=0; i<mtb.bs(index-1).numAtoms(); i++){
      AtomTopology at=mtb.bs(index-1).atom(i);
      at.setMass(gff.findMass(int(at.mass())));
      st.addAtom(at);
    }
    
    // polarisation
    for (int i=0; i<st.numAtoms(); i++) {
      st.atom(i).setPolarisable(mtb.bs(index-1).atom(i).isPolarisable());
      st.atom(i).setPolarisability(mtb.bs(index-1).atom(i).polarisability());
      st.atom(i).setCosCharge(mtb.bs(index-1).atom(i).cosCharge());
      st.atom(i).setDampingLevel(mtb.bs(index-1).atom(i).dampingLevel());
      st.atom(i).setDampingPower(mtb.bs(index-1).atom(i).dampingPower());
      st.atom(i).setPoloffsiteGamma(mtb.bs(index-1).atom(i).poloffsiteGamma());
      st.atom(i).setPoloffsiteI(mtb.bs(index-1).atom(i).poloffsiteI());
      st.atom(i).setPoloffsiteJ(mtb.bs(index-1).atom(i).poloffsiteJ());
    }
    
    ConstraintIterator cit(mtb.bs(index-1));
    for(;cit;++cit)
      st.addConstraint(cit());
    st.setSolvName(mtb.bs(index-1).solvName());

    sys.addSolvent(Solvent(st));
    
    // we have to determine still what is a H and what not
    for(int m=0; m<sys.numMolecules(); m++){
      sys.mol(m).topology().clearH();
      sys.mol(m).topology().setHmass(1.008);
    }

    // set the physical constants in the gff    
    gff.setFpepsi(mtb.Fpepsi());
    gff.setHbar(mtb.Hbar());
    gff.setSpdl(mtb.Spdl());
    gff.setBoltz(mtb.Boltz());

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
    title << "MAKE_TOP topology, using:" << endl;
    iter=args.lower_bound("build");
    to=args.upper_bound("build");
    for( ; iter!=to ; ++iter)
      title << iter->second << endl;
    iter=args.lower_bound("param");
    to=args.upper_bound("param");
    for( ; iter!=to ; ++iter)
      title << iter->second << endl;
 
    if(gff.ForceField()!="_no_FORCEFIELD_block_given_")
      title << endl << "Force-field code: "+gff.ForceField();

	ot.setTitle(title.str());

	ot.write(sys,gff);
    
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}
