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
#include "SimplePairlist.h"

#include <string>
#include <cassert>

#include "AtomSpecifier.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/Exclusion.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/Box.h"
#include "../bound/Boundary.h"
#include "../gmath/Vec.h"
#include "../gromos/Exception.h"

namespace utils{

  SimplePairlist::SimplePairlist(gcore::System &sys, 
				 bound::Boundary &pbc, 
				 double c){
    setSystem(sys);
    d_pbc = &pbc;
    d_cut2 = c*c;
    d_chargeGroupBased=false;
  }
  
  void SimplePairlist::setPbc(bound::Boundary &pbc)
  {
    d_pbc = &pbc;
  }
  
  void SimplePairlist::setCutOff(double c)
  {
    d_cut2 = c*c;
  }
  
  void SimplePairlist::setType(std::string s)
  {
    if(s=="ATOMIC")
      d_chargeGroupBased = false;
    else if(s=="CHARGEGROUP")
      d_chargeGroupBased = true;
    else
      throw gromos::Exception("SimplePairlist", 
			      "Pairlist type " + s + " unknown");
  }

  void SimplePairlist::setAtom(int m, int a)
  {
    if(m==-3)
      throw gromos::Exception("SimplePairlist", "you cannot use the function setAtom(m,a) for a virtual atom. Ask for help");
    if(m<0)
      d_atom = new SolventSpecAtom(*sys(), m, a);
    else
      d_atom = new SpecAtom(*sys(), m, a);
  }

  void SimplePairlist::setAtom(SpecAtom &s)
  {
    d_atom = &s;
  }
  
  void SimplePairlist::calc()
  {
    if(d_chargeGroupBased)
      calcCgb();
    else
      calcAtomic();
  }
  
  void SimplePairlist::calc(const AtomSpecifier &B, double cutmin)
  {
    if(d_chargeGroupBased)
      calcCgb(B, cutmin);
    else
      calcAtomic(B, cutmin);
  }

  void SimplePairlist::calcCgb()
  {
    gmath::Vec atom_i;
    if(d_atom->type() == spec_solvent && d_atom->atom() >= sys()->sol(0).numPos())
      throw gromos::Exception("SimplePairlist",
				"Not enough solvent atoms in system");

    atom_i = chargeGroupPosition(*d_atom);
    
    // gather the system to get the charge groups connected
    // d_pbc->gather();

    // now loop over all charge gropus and add those atoms that belong to
    // a charge group that is within d_cut2
    double d2;
    gmath::Vec v;
    
    //first solute
    for(int m=0; m<sys()->numMolecules(); m++){
      int a1=-1;
      int a2=0;
      
      
      while(a2<sys()->mol(m).numAtoms()){
	// search for the first chargeGroup==1
	for(; sys()->mol(m).topology().atom(a2).chargeGroup()!=1; a2++);
	v=chargeGroupPosition(m,a2);
	v=d_pbc->nearestImage(atom_i, v, sys()->box());
	d2=(atom_i - v).abs2();
	if(d2 <=d_cut2)
	  // now add everything from a1+1 to a2;
	  for(int i=a1; i<a2; i++) addAtom(m,i+1);
	a1 = a2;
	a2++;
      }
    }

    //and solvent
    int nsa = sys()->sol(0).topology().numAtoms();
    for(int i=0; i< sys()->sol(0).numPos(); i+=nsa){
      v = d_pbc->nearestImage(atom_i, sys()->sol(0).pos(i), sys()->box());
      d2 = (atom_i - v).abs2();
      if(d2 <= d_cut2)
	for(int j=0; j<nsa; j++) addAtom(-1,i+j);
    }
    
    // and remove the atom itself
    removeAtom(d_atom->mol(), d_atom->atom());
  }
  
  void SimplePairlist::calcAtomic()
  {
    
    if(d_atom->type() == spec_solvent &&
       d_atom->atom() >= sys()->sol(0).numPos())
      throw gromos::Exception("SimplePairlist",
				"Not enough solvent atoms in system");

    gmath::Vec atom_i = d_atom->pos();
    
    // now loop over all atoms and add those atoms that are within d_cut2
    double d2;
    gmath::Vec v;
    
    // first solute
    for(int m=0; m<sys()->numMolecules(); m++){
      for(int a=0; a<sys()->mol(m).numPos(); a++){
	v=d_pbc->nearestImage(atom_i, sys()->mol(m).pos(a), sys()->box());
	d2=(atom_i - v).abs2();
	if(d2<=d_cut2)
	  addAtom(m,a);
      }
    }
    // and solvent
    for(int i=0; i<sys()->sol(0).numPos(); i++){
      v=d_pbc->nearestImage(atom_i, sys()->sol(0).pos(i), sys()->box());
      d2=(atom_i - v).abs2();
      if(d2<=d_cut2)
	addAtom(-1,i);
    }
    // now remove the atom itself
    removeAtom(d_atom->mol(), d_atom->atom());
  }
  
  void SimplePairlist::calcAtomic(const AtomSpecifier &B, double cutmin)
  {
    // the position of the center atom
    gmath::Vec pos_i = d_atom->pos();
    
    // now loop over all atoms of group B and add those atoms that are within d_cut2
    double d2;
    double d2min = cutmin * cutmin;
    gmath::Vec v;
    for (unsigned int b = 0; b < B.size(); b++) {
      v = d_pbc->nearestImage(pos_i, B.pos(b), sys()->box());
      d2 = (pos_i - v).abs2();
      if (d2 < d_cut2 && d2 >= d2min) {
        addAtom(B.mol(b), B.atom(b));
      }
    }
    
    // now remove the atom itself
    removeAtom(d_atom->mol(), d_atom->atom());
  }
  
  void SimplePairlist::calcCgb(const AtomSpecifier &B, double cutmin)
  {
    // the position of the center charge group
    gmath::Vec pos_i;
    pos_i = chargeGroupPosition(*d_atom);

    // loop over all atoms of B and add those charge groups which are within d_cut2
    // it assumes that B contains complete charge groups only which must be tested
    // before calling this function!!!
    double d2;
    double d2min = cutmin * cutmin;
    gmath::Vec v;
    for (unsigned int b = 0; b < B.size(); b++) {
      int molNumB = B.mol(b);
      int atomNumB = B.atom(b);
      v = chargeGroupPosition(molNumB, atomNumB);
      v = d_pbc->nearestImage(pos_i, v, sys()->box());
      d2= (pos_i - v).abs2();
      // add the whole charge group B in case the distance of the two charge
      // groups is < sqrt(d_cut2)
      if(d2 < d_cut2 && d2 >= d2min) {
        // find the last atom of this charge group
        int lastAtom = atomNumB;
        for(;sys()->mol(molNumB).topology().atom(lastAtom).chargeGroup() != 1; lastAtom++);
        for(int a = atomNumB; a <= lastAtom; a++) {
          addAtom(molNumB, a);
          b++;
        }
        b--;
      }
    }
    // and remove the atom itself
    removeAtom(d_atom->mol(), d_atom->atom());
  }
  
  gmath::Vec SimplePairlist::chargeGroupPosition(int m, int a)
  {
    gmath::Vec v(0.0,0.0,0.0);
    // solvent
    if(m<0){
      int i=a/sys()->sol(0).topology().numAtoms();
      i*=sys()->sol(0).topology().numAtoms();
      return sys()->sol(0).pos(i);
    }
    int begin=a-1, end=a;
    if(a>0)
      for(begin=a-1;
	  begin>=0 && sys()->mol(m).topology().atom(begin).chargeGroup()!=1; 
	  begin--);
    
    for(end=a;
	sys()->mol(m).topology().atom(end).chargeGroup()!=1;
	end++);
    
    // charge group goes from begin+1 to end
    for(int k=begin+1; k<=end; k++)
      v += d_pbc->nearestImage(sys()->mol(m).pos(begin+1),
			       sys()->mol(m).pos(k),
			       sys()->box());
    return v/(end-begin);
  }
  
  gmath::Vec SimplePairlist::chargeGroupPosition(SpecAtom &s)
  {
    int m=s.mol();
    int a=s.atom();
    
    gmath::Vec v(0.0,0.0,0.0);
    // solvent
    if(s.type() == spec_solvent){
      int i=a/sys()->sol(0).topology().numAtoms();
      i*=sys()->sol(0).topology().numAtoms();
      return sys()->sol(0).pos(i);
    }
    // virtual atom
    if(s.type() == spec_virtual){
      return s.pos();
    }
    
    int begin=a-1, end=a;
    if(a>0)
      for(begin=a-1;
	  begin>=0 && sys()->mol(m).topology().atom(begin).chargeGroup()!=1; 
	  begin--);
    
    for(end=a;
	sys()->mol(m).topology().atom(end).chargeGroup()!=1;
	end++);
    
    // charge group goes from begin+1 to end
    for(int k=begin+1; k<=end; k++)
      v += d_pbc->nearestImage(sys()->mol(m).pos(begin+1),
			       sys()->mol(m).pos(k),
			       sys()->box());
    return v/(end-begin);
  }

  void SimplePairlist::removeExclusions()
  {
    // of course it is not effective to first add and later remove
    // the excluded atoms, but if you only want a list of atoms within a
    // cutoff then you just do not call this function

    // check whether we are looking at a solvent
    if(d_atom->type() == spec_solvent){
      int nsa=sys()->sol(0).topology().numAtoms();
      int first=d_atom->atom()/nsa;
      first *= nsa;
      for(int i=0; i<nsa; i++) removeAtom(-1, first+i);
    }
    else{
      
    // loop over all solute atoms before d_atom
      for(int a=0; a<d_atom->atom(); a++)
	for(int i=0; 
	    i< sys()->mol(d_atom->mol()).topology().atom(a).exclusion().size();
	    i++)
	  if(d_atom->atom() == 
	     sys()->mol(d_atom->mol()).topology().atom(a).exclusion().atom(i))
	    removeAtom(d_atom->mol(),a);
      // and remove all excluded atoms of d_a
      for(int i=0; 
	  i < sys()->mol(d_atom->mol()).topology().atom(d_atom->atom()).exclusion().size();
	  i++)
	removeAtom(d_atom->mol(), 
	   sys()->mol(d_atom->mol()).topology().atom(d_atom->atom()).exclusion().atom(i));
    }
  }

  void SimplePairlist::remove14Exclusions()
  {
    // of course it is not effective to first add and later remove
    // the excluded atoms, but if you only want a list of atoms within a
    // cutoff then you just do not call this function

    // check whether we are looking at a solvent
    if(d_atom->type() == spec_solvent){
      return;
    }
    else{
      
    // loop over all solute atoms before d_atom
      for(int a=0; a<d_atom->atom(); a++)
	for(int i=0; 
	    i< sys()->mol(d_atom->mol()).topology().atom(a).exclusion14().size();
	    i++)
	  if(d_atom->atom() == 
	     sys()->mol(d_atom->mol()).topology().atom(a).exclusion14().atom(i))
	    removeAtom(d_atom->mol(),a);
      // and remove all excluded atoms of d_a
      for(int i=0; 
	  i < sys()->mol(d_atom->mol()).topology().atom(d_atom->atom()).exclusion14().size();
	  i++)
	removeAtom(d_atom->mol(), 
	   sys()->mol(d_atom->mol()).topology().atom(d_atom->atom()).exclusion14().atom(i));
    }
  }

}



     
