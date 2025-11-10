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

// gcore_LinearTopology.cc

#include "LinearTopology.h"

#include <algorithm>
#include <cassert>
#include <string>
#include <set>
#include <utility>
#include <vector>
#include <map>

#include "Molecule.h"
#include "VirtualAtoms.h"
#include "AtomTopology.h"
#include "Exclusion.h"
#include "Bond.h"
#include "Constraint.h"
#include "Angle.h"
#include "Improper.h"
#include "Dihedral.h"
#include "CrossDihedral.h"
#include "LJException.h"
#include "System.h"
#include "MoleculeTopology.h"
#include "LJException.h"
#include "../gromos/Exception.h"

using namespace std;
using namespace gcore;

LinearTopology::LinearTopology(){};

LinearTopology::LinearTopology(gcore::System &sys)
{
  int lastAtom=0;
  int lastResidue=0;

  for(int m=0; m< sys.numMolecules(); m++){
    for(int a=0; a<sys.mol(m).numAtoms(); a++){
      AtomTopology at=sys.mol(m).topology().atom(a);
      Exclusion ex, ex14;
      for(int e=0; e< sys.mol(m).topology().atom(a).exclusion().size(); e++){
	ex.insert(sys.mol(m).topology().atom(a).
		  exclusion().atom(e)+lastAtom);
      }
      for(int e=0; e< sys.mol(m).topology().atom(a).exclusion14().size(); 
	  e++){
	ex14.insert(sys.mol(m).topology().atom(a).
		  exclusion14().atom(e)+lastAtom);
      }
      at.setExclusion(ex);
      at.setExclusion14(ex14);
      
      d_atom.push_back(at);
      d_resmap[a+lastAtom]=sys.mol(m).topology().resNum(a)+lastResidue;
    }
    for(int i=0; i<sys.mol(m).topology().numRes(); i++)
      d_resname.push_back(sys.mol(m).topology().resName(i));

    BondIterator bi(sys.mol(m).topology());
    for(; bi; ++bi){
      Bond b=bi();
      b[0]+=lastAtom; b[1]+=lastAtom;
      d_bond.insert(b);
    }

    set<Constraint>::const_iterator iter = sys.mol(m).topology().constraints().begin(), 
                                     to = sys.mol(m).topology().constraints().end();
    for(; iter != to; ++iter){
      Constraint b=*iter;
        b[0]+=lastAtom; b[1]+=lastAtom;
        d_constraints.insert(b);
    }

    BondDipoleIterator bdi(sys.mol(m).topology());
    for(; bdi; ++bdi){
      Bond b=bdi();
      b[0]+=lastAtom; b[1]+=lastAtom;
      d_dipole_bond.insert(b);
    }    

    AngleIterator ai(sys.mol(m).topology());
    for(; ai; ++ai){
      Angle a=ai();
      a[0]+=lastAtom; a[1]+=lastAtom; a[2]+=lastAtom;
      d_angle.insert(a);
    }

    DihedralIterator di(sys.mol(m).topology());
    for(; di; ++di){
      Dihedral d=di();
      d[0]+=lastAtom; d[1]+=lastAtom; d[2]+=lastAtom; d[3]+=lastAtom;
      d_dihedral.insert(d);
    }

    CrossDihedralIterator cdi(sys.mol(m).topology());
    for (; cdi; ++cdi) {
      CrossDihedral cd = cdi();
      for (unsigned int i = 0; i < 8; ++i)
        cd[i] += lastAtom;
      d_crossdihedral.insert(cd);
    }

    ImproperIterator ii(sys.mol(m).topology());
    for(; ii; ++ii){
      Improper i=ii();
      i[0]+=lastAtom; i[1]+=lastAtom; i[2]+=lastAtom; i[3]+=lastAtom;
      d_improper.insert(i);
    }

    LJExceptionIterator lji(sys.mol(m).topology());
    for(; lji; ++lji){
      LJException lj=lji();
      lj[0]+=lastAtom;
      lj[1]+=lastAtom;
      d_ljexception.insert(lj);
    }
    
    lastResidue += sys.mol(m).topology().numRes();
    lastAtom    += sys.mol(m).numAtoms();
  }
}

LinearTopology::~LinearTopology(){};
gcore::System LinearTopology::parse()
{
  gcore::System sys;
  parse(sys);
  return sys;
}

void LinearTopology::parse(gcore::System &sys)
{
  // largely copied from InTopology
  // res_corr is a correction to the residue number for
  // weird cases if a molecule is split in to two, in the middle of
  // a residue
  unsigned int atomCounter=0;
  set<Bond>::const_iterator bi=d_bond.begin();
  set<Constraint>::const_iterator ci=d_constraints.begin();
  set<Bond>::const_iterator bdi=d_dipole_bond.begin();
  set<Angle>::const_iterator ai=d_angle.begin();
  set<Improper>::const_iterator ii=d_improper.begin();
  set<Dihedral>::const_iterator di=d_dihedral.begin();
  set<CrossDihedral>::const_iterator cdi=d_crossdihedral.begin();
  set<LJException>::const_iterator lji=d_ljexception.begin();
  
//  int prevMol=0, lastAtom=0, prevMolRes=0, resCorr=0;
  int prevMol=0, lastAtom=0, prevMolRes=0, resCorr=0;
  
  MoleculeTopology *mt;

  // std::cerr << "parsing the linear topology!" << std::endl;
  // std::cerr << "d_atom size = " << d_atom.size() << std::endl;
  
  while(atomCounter < d_atom.size()){

    // std::cerr << "a new molecule" << std::endl;
    
    mt=new MoleculeTopology();

    // Detect the last atom of the first molecule & add bonds:
    // std::cerr << "detect last atom" << std::endl;
    
    for( ;bi != d_bond.end() && (*bi)[0] <= lastAtom; ++bi){
      Bond bond = *bi;
      if(bond[1]>lastAtom){
	lastAtom=bond[1];
	if (unsigned(lastAtom) >= d_atom.size())
	  throw gromos::Exception("LinearTopology", "bonds between non-existing atoms");
      }
      bond[0] -= prevMol; bond[1] -= prevMol;
      if(bond[0]>=0 && bond[1]>=0)
        mt->addBond(bond);
    }
    lastAtom++;
    // std::cerr << "last atom = " << lastAtom << std::endl;

    // add DipoleBonds
    /*
     * This is a check to add the atoms that do not have a GROMOS bond
     * but have dipole bonds
     * 
     */
    bool lastAtomChanged = false;
    for( ;bdi != d_dipole_bond.end() && (*bdi)[0] < lastAtom; ++bdi)
    {
      Bond bond = *bdi;
      if(bond[1]>=lastAtom){
        lastAtom=bond[1];
        lastAtomChanged=true;
	if (unsigned(lastAtom) >= d_atom.size())
	  throw gromos::Exception("LinearTopology", "dipole bonds between non-existing atoms");
      }
      bond[0] -= prevMol; bond[1] -= prevMol;
      if(bond[0]>=0 && bond[1]>=0)
        mt->addDipoleBond(bond);
    }
    if(lastAtomChanged)
        lastAtom++;

    // add Atoms
    //std::cerr<< "sysres "<< sys.mol(0).topology().numRes();
    for(; int(atomCounter) < lastAtom; atomCounter++){
       
      //std::cerr << "adding atom " << atomCounter << std::endl;
      // convert poloffsite atom numbers to be relative to the molecule
      if (d_atom[atomCounter].isPolarisable()) {
        int poloffsite_i=d_atom[atomCounter].poloffsiteI();
        int poloffsite_j=d_atom[atomCounter].poloffsiteJ();
        d_atom[atomCounter].setPoloffsiteI(poloffsite_i-prevMol);
        d_atom[atomCounter].setPoloffsiteJ(poloffsite_j-prevMol);
      }

      mt->addAtom(d_atom[atomCounter]);
      
      // adapt exclusions:
      Exclusion *e;
      e=new Exclusion();
      for (int l=0;l<d_atom[atomCounter].exclusion().size();++l){
        e->insert(d_atom[atomCounter].exclusion().atom(l) - prevMol);
      }
      mt->atom(mt->numAtoms()-1).setExclusion(*e);
      delete e;
      e=new Exclusion();
      for (int l=0;l<d_atom[atomCounter].exclusion14().size();++l){
        e->insert(d_atom[atomCounter].exclusion14().atom(l)-prevMol);
      }
      mt->atom(mt->numAtoms()-1).setExclusion14(*e);
      delete e;
      
      int resn=d_resmap[atomCounter]-prevMolRes;

       if(resn+resCorr<0){
           resCorr -= resn;
       };      
      mt->setResNum(atomCounter-prevMol,resn+resCorr);
      mt->setResName(resn+resCorr,d_resname[resn+prevMolRes]);
    }
    prevMolRes+=mt->numRes();

    // add DipoleBonds
//    for( ;bdi != d_dipole_bond.end() && (*bdi)[0] < lastAtom; ++bdi)
//    {
//      Bond bond = *bdi;
//      if(bond[1]>lastAtom){
//	lastAtom=bond[1];
//	if (unsigned(lastAtom) >= d_atom.size())
//	  throw gromos::Exception("LinearTopology", "dipole bonds between non-existing atoms");
//      }
//      bond[0] -= prevMol; bond[1] -= prevMol;
//      if(bond[0]>=0 && bond[1]>=0)
//        mt->addDipoleBond(bond);
//    }
    
    // add Constraints
    for( ; ci != d_constraints.end() && (*ci)[0] < lastAtom; ++ci){
      Constraint constraint = *ci;
      constraint[0] -= prevMol; constraint[1] -= prevMol;
      if(constraint[0]>=0 && constraint[1]>=0)
        mt->addConstraint(constraint);
    }
    
    // add Angles
    for( ; ai != d_angle.end() && (*ai)[0] < lastAtom; ++ai){
      Angle angle = *ai;
      angle[0] -= prevMol; angle[1] -= prevMol; angle[2] -= prevMol;
      if(angle[0]>=0 && angle[1]>=0 && angle[3]>=0)
        mt->addAngle(angle);
    }
    
    // add Dihedrals
    for( ; di != d_dihedral.end() && (*di)[0] < lastAtom; ++di)
    {
      Dihedral dihedral = *di;
      dihedral[0] -= prevMol; dihedral[1] -= prevMol;
      dihedral[2] -= prevMol; dihedral[3] -= prevMol;
      if(dihedral[0]>=0 && dihedral[1]>=0 && dihedral[2]>=0 && dihedral[3]>=0)
        mt->addDihedral(dihedral);
    }

    // add CrossDihedrals
    for( ; cdi != d_crossdihedral.end() && (*cdi)[0] < lastAtom; ++cdi)
    {
      CrossDihedral crossdihedral = *cdi;
      bool positive = true;
      for(unsigned int i = 0; i < 8; ++i) {
        crossdihedral[i] -= prevMol;
        positive = positive && (crossdihedral[i] > 0 ? true : false);
      }
      if(positive)
        mt->addCrossDihedral(crossdihedral);
    }
    
    // add Impropers
    for( ; ii != d_improper.end() && (*ii)[0] < lastAtom; ++ii){
      Improper improper = *ii;
      improper[0] -= prevMol; improper[1] -= prevMol;
      improper[2] -= prevMol; improper[3] -= prevMol;
      if(improper[0]>=0 && improper[1]>=0 && improper[2]>=0 && improper[3]>=0)
        mt->addImproper(improper); 
    }

    // add LJ exceptions
    for( ; lji != d_ljexception.end() && (*lji)[0] < lastAtom; ++lji){
      LJException ljexception = *lji;
      ljexception[0] -= prevMol;
      ljexception[1] -= prevMol;
      if(ljexception[0]>=0 && ljexception[1]>=0)
        mt->addLJException(ljexception);
    }

    // add the molecule to the system.
    sys.addMolecule(Molecule(*mt));
    delete mt;
    prevMol=lastAtom;
  }
  // std::cerr << "LinearTopology::parse done" << std::endl;
}


void LinearTopology::get14s()
{
  int na=d_atom.size();

  for(int i=0; i<na; i++){
    set<int> first, second, third;
    set<Bond>::const_iterator bi1 = d_bond.begin(), bi2, bi3;
    for(bi1=d_bond.begin(); bi1 != d_bond.end(); ++bi1){
      if(i == (*bi1)[0]) first.insert((*bi1)[1]);
      if(i == (*bi1)[1]) first.insert((*bi1)[0]);
    }
    for(set<int>::const_iterator iter=first.begin(), to=first.end();
        iter!=to; ++iter){
      for(bi2=d_bond.begin(); bi2 != d_bond.end(); ++bi2){
        if(*iter == (*bi2)[0]) second.insert((*bi2)[1]);
        if(*iter == (*bi2)[1]) second.insert((*bi2)[0]);
      }
    }
    for(set<int>::const_iterator iter=second.begin(), to=second.end();
        iter!=to; ++iter){ 
      for(bi3=d_bond.begin(); bi3 != d_bond.end(); ++bi3){
        if(*iter == (*bi3)[0]) third.insert((*bi3)[1]);
        if(*iter == (*bi3)[1]) third.insert((*bi3)[0]);
      }
    }
    Exclusion e;
    for(set<int>::const_iterator iter=third.begin(), to=third.end();
        iter!=to; ++iter){
      if(i<*iter&&
         !first.count(*iter)&&
         !second.count(*iter)){
        int excl=0;
        for(int k=0; k < d_atom[i].exclusion().size();k++)
          if(*iter == d_atom[i].exclusion().atom(k)) excl=1;
        if(!excl) e.insert(*iter);
      }
    }
    d_atom[i].setExclusion14(e);
  }
}

std::vector<int> LinearTopology::removeAtoms()
{
  set<int> rem;
  vector<int> ren;
  int corr=0;
  for(unsigned int i=0; i< d_atom.size();i++){

    if ( d_atom[i].iac() < 0 ) {
      rem.insert(i);
      corr++;
      ren.push_back(-1);
    }
    else{
      ren.push_back(i-corr);
    }
  }
  if ( rem.size() == 0 ) return ren;
  // add four more to ren, in order to have a buffer
  // and why did we need this? I think for cyclization in maketop
  for(int i=0; i<6; i++)
    ren.push_back(d_atom.size()+i-corr);
  
  // process the properties one at a time 
  _reduceResidues(rem, ren);  //has to be done before reduce atoms, as otherwise there are less atoms,than in the resMap. bschroed
  _reduceAtoms(rem, ren);
  _reduceBonds(rem, ren);
  _reduceConstraints(rem, ren);
  _reduceDipoleBonds(rem, ren);
  _reduceAngles(rem, ren);
  _reduceImpropers(rem, ren);
  _reduceDihedrals(rem, ren);
  _reduceCrossDihedrals(rem, ren);
  _reduceLJExceptions(rem,ren);

  return ren;
}


void LinearTopology::_reduceAtoms(std::set<int> &rem, std::vector<int> &ren)
{
  int count=0;
  for(vector<AtomTopology>::iterator iter=d_atom.begin();
      iter!=d_atom.end();){
    if (rem.count(count)){
        d_atom.erase(iter);
    }
    else{
      Exclusion e;
      for(int j=0; j < iter->exclusion().size(); j++){
        if(!rem.count(iter->exclusion().atom(j)))
          e.insert(ren[iter->exclusion().atom(j)]);
      }
      iter->setExclusion(e);
      ++iter;
    }
    count++;
  }
}

void LinearTopology::_reduceResidues(std::set<int> &rem, std::vector<int> &ren)
{
    /*
     * This function is renumbering the residue-atom linkeage in red_top. 
     * adapted  bschroed
     */
  
  //DEBUG BEFORE - remove maybe? bschroed
  /*
  std::cerr << "BEFORE: ResiNamesize "<< d_resname.size() << std::endl;//todo: remove!@
  std::cerr << "BEFORE: Total ATOMS "<< d_atom.size()<< std::endl; //todo: remove!@
  std::cerr << "BEFORE: ResMap size (should be equal to total atomsize!) "<< d_resmap.size()<< "\n\n"; //todo: remove!@
  */

  int lastRes=-1;   //old_residue number
  int resNum=-1;    //new residue num
  map<int, int> tempMap;
  vector<string> tempNames;
  
  for(int i=0; i < d_atom.size(); i++){
    if(!rem.count(i)){  //if the atom is not in the remove list.
      if(d_resmap[i] != lastRes){   //a new residue?
        lastRes=d_resmap[i];
        tempNames.push_back(d_resname[lastRes]);
        resNum++;       
      }
      tempMap.insert(std::pair<int, int>(ren[i], resNum));  //add entry
    }
   }
  d_resmap  = tempMap;
  d_resname = tempNames;
  
  //DEBUG AFTER - remove maybe? bschroed
  /*
  for (auto x: tempMap){
      std::cerr << "TEMP MAP\t" << x.first << "/"<< x.second << std::endl;
  }

  for (auto x: tempNames){
      std::cerr << "TEMP NAMES\t" << x<< std::endl;
  }

  std::cerr << "AFTER: Total Res Num: " << d_resmap.size()<< std::endl;//todo: remove!@
  std::cerr << "AFTER: reduce ResiNamesize "<< d_resname.size()<< std::endl; //todo: remove!@
  */
}

void LinearTopology::_reduceBonds(std::set<int> &rem, std::vector<int> &ren)
{
  //these are a set. Changing them while looping over them will change the
  // order during the loop. Rather create a new set and copy over...
  set<Bond> newBonds;
  set<Bond>::const_iterator iter = d_bond.begin(), to=d_bond.end();
  for(; iter != to; ++iter){
    if(rem.count((*iter)[0]) == 0 && rem.count((*iter)[1]) == 0){
      Bond b(ren[(*iter)[0]], ren[(*iter)[1]]);
      b.setType(iter->type());
      newBonds.insert(b);
    }
  }
  d_bond = newBonds;
}

void LinearTopology::_reduceConstraints(std::set<int> &rem, std::vector<int> &ren)
{
  //these are a set. Changing them while looping over them will change the
  // order during the loop. Rather create a new set and copy over...
  set<Constraint> newConstraints;
  set<Constraint>::const_iterator iter = d_constraints.begin(), to=d_constraints.end();
  for(; iter != to; ++iter){
    if(rem.count((*iter)[0]) == 0 && rem.count((*iter)[1]) == 0){
      Constraint b(ren[(*iter)[0]], ren[(*iter)[1]]);
      b.setType(iter->bondtype());
      b.setDist(iter->dist());
      newConstraints.insert(b);
    }
  }
  d_constraints = newConstraints;
}

void LinearTopology::_reduceDipoleBonds(std::set<int> &rem, std::vector<int> &ren)
{
  //these are a set. Changing them while looping over them will change the
  // order during the loop. Rather create a new set and copy over...
  set<Bond> newBonds;
  set<Bond>::const_iterator iter = d_dipole_bond.begin(), to=d_dipole_bond.end();
  for(; iter != to; ++iter){
    if(rem.count((*iter)[0]) == 0 && rem.count((*iter)[1]) == 0){
      Bond b(ren[(*iter)[0]], ren[(*iter)[1]]);
      b.setType(iter->type());
      newBonds.insert(b);
    }
  }
  d_dipole_bond = newBonds;
}

void LinearTopology::_reduceAngles(std::set<int> &rem, std::vector<int> &ren)
{
  //these are a set. Changing them while looping over them will change the
  // order during the loop. Rather create a new set and copy over...
  set<Angle> newAngles;
  set<Angle>::const_iterator iter = d_angle.begin(), to=d_angle.end();
  for(; iter != to; ++iter){
    if(rem.count((*iter)[0]) == 0 && rem.count((*iter)[1]) == 0 &&
       rem.count((*iter)[2]) == 0){
      Angle a(ren[(*iter)[0]], ren[(*iter)[1]], ren[(*iter)[2]]);
      a.setType(iter->type());
      newAngles.insert(a);
    }
  }
  d_angle = newAngles;
}

void LinearTopology::_reduceImpropers(std::set<int> &rem, std::vector<int> &ren)
{
  //these are a set. Changing them while looping over them will change the
  // order during the loop. Rather create a new set and copy over...
  set<Improper> newImpropers;
  set<Improper>::const_iterator iter = d_improper.begin(), to=d_improper.end();
  for(; iter != to; ++iter){
    if(rem.count((*iter)[0]) == 0 && rem.count((*iter)[1]) == 0 &&
       rem.count((*iter)[2]) == 0 && rem.count((*iter)[3]) == 0){
      Improper i(ren[(*iter)[0]], ren[(*iter)[1]], 
		 ren[(*iter)[2]], ren[(*iter)[3]]);
      i.setType(iter->type());
      newImpropers.insert(i);
    }
  }
  d_improper = newImpropers;
}


void LinearTopology::_reduceDihedrals(std::set<int> &rem, std::vector<int> &ren)
{
   //these are a set. Changing them while looping over them will change the
  // order during the loop. Rather create a new set and copy over...
  set<Dihedral> newDihedrals;
  set<Dihedral>::const_iterator iter = d_dihedral.begin(), to=d_dihedral.end();
  for(; iter != to; ++iter){
    if(rem.count((*iter)[0]) == 0 && rem.count((*iter)[1]) == 0 &&
       rem.count((*iter)[2]) == 0 && rem.count((*iter)[3]) == 0){
      Dihedral i(ren[(*iter)[0]], ren[(*iter)[1]], 
		 ren[(*iter)[2]], ren[(*iter)[3]]);
      i.setType(iter->type());
      newDihedrals.insert(i);
    }
  }
  d_dihedral = newDihedrals;
}

void LinearTopology::_reduceCrossDihedrals(std::set<int> &rem, std::vector<int> &ren)
{
   //these are a set. Changing them while looping over them will change the
  // order during the loop. Rather create a new set and copy over...
  set<CrossDihedral> newCrossDihedrals;
  set<CrossDihedral>::const_iterator iter = d_crossdihedral.begin(), to=d_crossdihedral.end();
  for(; iter != to; ++iter){
    if(rem.count((*iter)[0]) == 0 && rem.count((*iter)[1]) == 0 &&
       rem.count((*iter)[2]) == 0 && rem.count((*iter)[3]) == 0 &&
       rem.count((*iter)[4]) == 0 && rem.count((*iter)[5]) == 0 &&
       rem.count((*iter)[6]) == 0 && rem.count((*iter)[7]) == 0){
      CrossDihedral i(ren[(*iter)[0]], ren[(*iter)[1]],
		 ren[(*iter)[2]], ren[(*iter)[3]],
                 ren[(*iter)[4]], ren[(*iter)[5]],
                 ren[(*iter)[6]], ren[(*iter)[7]]);
      i.setType(iter->type());
      newCrossDihedrals.insert(i);
    }
  }
  d_crossdihedral = newCrossDihedrals;
}

void LinearTopology::_reduceLJExceptions(std::set<int> &rem, std::vector<int> &ren)
{
  //these are a set. Changing them while looping over them will change the
  // order during the loop. Rather create a new set and copy over...
  set<LJException> newLJExceptions;
  set<LJException>::const_iterator iter = d_ljexception.begin(), to=d_ljexception.end();
  for(; iter != to; ++iter){
    if(rem.count((*iter)[0]) == 0 && rem.count((*iter)[1]) == 0){
      LJException lj(ren[(*iter)[0]],ren[(*iter)[1]]);
      lj.setType(iter->type());
      newLJExceptions.insert(lj);
    }
  }
  d_ljexception = newLJExceptions;
}

void LinearTopology::moveAtoms(std::vector<std::pair<int, int> > moveatoms) {

  // loop over from/to pairs
  std::vector<std::pair<int, int> >::iterator p_it = moveatoms.begin(),
                                              p_to = moveatoms.end();
  for (;p_it != p_to;p_it++) {
    unsigned int move_from=p_it->first-1;
    unsigned int move_to=p_it->second-1;
    int upper=std::max(move_from,move_to);
    
    if (d_resmap[move_from] != d_resmap[move_to])  {
        // restrict to moving atoms within one residue, otherwise
        // it is not defined which residue they should be assigned to
        throw gromos::Exception("LinearTopology::moveAtoms", "both atoms of a pair have to be from the same residue!");
    } else if (move_from > d_atom.size() || move_to > d_atom.size()) { 
        throw gromos::Exception("LinearTopology::moveAtoms", "atom number out of range!");
    }
    
    vector<AtomTopology>::iterator firstatom = d_atom.begin();
    AtomTopology at(d_atom[move_from]);
    
    std::map<int,int> change_map;
    // insert/remove atoms and fill map that connects old and new indices
    change_map[move_from]=move_to;
    if (move_from > move_to) {
        for (unsigned int ii=move_to; ii < move_from; ii++) {
          change_map[ii]=ii+1;
        }
        d_atom.erase(firstatom+move_from);
        d_atom.insert(firstatom+move_to, at);
    } else if (move_from < move_to) {
        for (unsigned int ii=move_from+1; ii <= move_to; ii++) {
          change_map[ii]=ii-1;
        }
        d_atom.insert(firstatom+move_to, at);
        d_atom.erase(firstatom+move_from);
    }

    //for (std::map<int, int>::iterator mit = change_map.begin(), mto=change_map.end(); mit != mto; mit++) {
    //  std::cerr<< mit->first << " - " << mit->second << std::endl;
    //}
    
    // renumber exclusions, 14excl
    for (vector<AtomTopology>::iterator atom_it = firstatom; atom_it <= firstatom+upper ; atom_it++) {
      Exclusion tempexcl, tempexcl14;
      for (int jj=0; jj < atom_it->exclusion().size(); jj++) { 
        int ex = atom_it->exclusion().atom(jj);
        if (change_map.count(ex)) { tempexcl.insert(change_map[ex]);
        }
        else tempexcl.insert(ex);
      }
      atom_it->exclusion()=tempexcl;
    
      for (int jj=0; jj < atom_it->exclusion14().size(); jj++) { 
        int ex = atom_it->exclusion14().atom(jj);
        if (change_map.count(ex))  tempexcl14.insert(change_map[ex]);
        else tempexcl14.insert(ex);
      }
      atom_it->exclusion14()=tempexcl14;
    }      
  
    // move exclusions to first atom involving them 
    for (int atom_i = 0; atom_i<=upper; atom_i++) {
      for (int jj=0; jj < d_atom[atom_i].exclusion().size(); jj++) { 
        int ex = d_atom[atom_i].exclusion().atom(jj);
        if (ex < atom_i) {
          d_atom[ex].exclusion().insert(atom_i);
          d_atom[atom_i].exclusion().erase(ex);
        }
      }
    }  
    for (int ii = 0; ii<=upper; ii++) {
      for (int jj=0; jj < d_atom[ii].exclusion14().size(); jj++) { 
        int ex = d_atom[ii].exclusion14().atom(jj);
        if (ex < ii) {
          d_atom[ex].exclusion14().insert(ii);
          d_atom[ii].exclusion14().erase(ex);
        }
      }
    }          
  
  
    // renumber bonds, angles, dihedrals    
    // I can not change existing bonds because they are const, so I have
    // to create a new vector
    set<Bond> newBonds;
    for(set<Bond>::iterator iter = d_bond.begin(), to=d_bond.end(); iter != to; ++iter){
      int a[2];
      for (int ii=0; ii < 2;ii++) { 
        a[ii]=(*iter)[ii];
        if (change_map.count((*iter)[ii])) {
          a[ii]=change_map[(*iter)[ii]];
        }
      }
      Bond b(a[0],a[1]);
      b.setType(iter->type());
      newBonds.insert(b);
    }
    d_bond = newBonds;
    newBonds.clear();
    for(set<Bond>::iterator iter = d_dipole_bond.begin(), to=d_dipole_bond.end(); iter != to; ++iter){
      int a[2];
      for (int ii=0; ii < 2;ii++) { 
        a[ii]=(*iter)[ii];
        if (change_map.count((*iter)[ii])) {
          a[ii]=change_map[(*iter)[ii]];
        }
      }
      Bond b(a[0],a[1]);
      b.setType(iter->type());
      newBonds.insert(b);
    }
    d_dipole_bond = newBonds;

    set<Constraint> newConstraints;
    for(set<Constraint>::iterator iter = d_constraints.begin(), to=d_constraints.end(); iter != to; ++iter){
      int a[2];
      for (int ii=0; ii < 2;ii++) { 
        a[ii]=(*iter)[ii];
        if (change_map.count((*iter)[ii])) {
          a[ii]=change_map[(*iter)[ii]];
        }
      }
      Constraint b(a[0],a[1]);
      b.setType(iter->bondtype());
      newConstraints.insert(b);
    }
    d_constraints = newConstraints;
    
    set<Angle> newAngles;  
    for(set<Angle>::const_iterator iter = d_angle.begin(), to=d_angle.end(); iter != to; ++iter){
      int a[3];
      for (int ii=0; ii < 3;ii++) { 
        a[ii]=(*iter)[ii];
        if (change_map.count((*iter)[ii])) {
          a[ii]=change_map[(*iter)[ii]];
        }
      }
      Angle b(a[0],a[1],a[2]);
      b.setType(iter->type());
      newAngles.insert(b);
    }
    d_angle = newAngles;
  
    set<Improper> newImpropers;  
    for(set<Improper>::const_iterator iter = d_improper.begin(), to=d_improper.end(); iter != to; ++iter){
      int a[4];
      for (int ii=0; ii < 4;ii++) { 
        a[ii]=(*iter)[ii];
        if (change_map.count((*iter)[ii])) {
          a[ii]=change_map[(*iter)[ii]];
        }
      }
      Improper b(a[0],a[1],a[2],a[3]);
      b.setType(iter->type());
      newImpropers.insert(b);
    }
    d_improper = newImpropers;
  
    set<Dihedral> newDihedrals;  
    for(set<Dihedral>::const_iterator iter = d_dihedral.begin(), to=d_dihedral.end(); iter != to; ++iter){
      int a[4];
      for (int ii=0; ii < 4;ii++) { 
        a[ii]=(*iter)[ii];
        if (change_map.count((*iter)[ii])) {
          a[ii]=change_map[(*iter)[ii]];
        }
      }
      Dihedral b(a[0],a[1],a[2],a[3]);
      b.setType(iter->type());
      newDihedrals.insert(b);
    }
    d_dihedral = newDihedrals;  
    set<CrossDihedral> newCrossDihedrals;  
    for(set<CrossDihedral>::const_iterator iter = d_crossdihedral.begin(), to=d_crossdihedral.end(); iter != to; ++iter){
      int a[8];
      for (int ii=0; ii < 8;ii++) { 
        a[ii]=(*iter)[ii];
        if (change_map.count((*iter)[ii])) {
          a[ii]=change_map[(*iter)[ii]];
        }
      }
      CrossDihedral b(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7]);
      b.setType(iter->type());
      newCrossDihedrals.insert(b);
    }
    d_crossdihedral = newCrossDihedrals;
  
    set<LJException> newLJExceptions;  
    for(set<LJException>::const_iterator iter = d_ljexception.begin(), to=d_ljexception.end(); iter != to; ++iter){
      int a[2];
      for (int ii=0; ii < 2;ii++) { 
        a[ii]=(*iter)[ii];
        if (change_map.count((*iter)[ii])) {
          a[ii]=change_map[(*iter)[ii]];
        }
      }
      LJException b(a[0],a[1]);
      b.setType(iter->type());
      newLJExceptions.insert(b);
    }
    d_ljexception = newLJExceptions;
  
  } // end loop over pairs
  
}
