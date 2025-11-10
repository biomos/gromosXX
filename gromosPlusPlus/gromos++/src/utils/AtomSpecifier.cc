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
#include "AtomSpecifier.h"

#include <cassert>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <string>
#include <vector>

#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gmath/Vec.h"

using namespace gcore;
using namespace std;

utils::AtomSpecifier::AtomSpecifier() {
  d_sys=NULL;
  d_nsm=-1;
  d_not_atoms = NULL;
}

std::string utils::SpecAtom::toString()const
{
  std::ostringstream os;
  os << d_mol+1 << ":" << d_atom+1;
  return os.str();
}

void utils::AtomSpecifier::_appendAtom(int m, int a)
{
  // check whether it is already in
  if(findAtom(m, a) == -1){
    if (m < 0)
      d_specatom.push_back(new SolventSpecAtom(*d_sys, m, a));
    else
      d_specatom.push_back(new SpecAtom(*d_sys, m, a));
  }
}

void utils::AtomSpecifier::_appendSolvent(int m, int a)const
{
  // check whether it is already in
  assert(m<0);
  if(findAtom(m, a) == -1){
    // if (m < 0)
    d_specatom.push_back(new SolventSpecAtom(*d_sys, m, a));
    // else
    // d_specatom.push_back(new SpecAtom(*d_sys, m, a));
  }
}

bool utils::AtomSpecifier::_checkName(int m, int a, std::string s)const
{
  std::string::size_type iterator=s.find('?');
  std::string name_in_topo;
  // take in the following three lines to allow for solvent as well
  if(m<0)
    name_in_topo=d_sys->sol(0).topology().atom(a).name().substr(0, iterator);
  else
    name_in_topo=d_sys->mol(m).topology().atom(a).name().substr(0, iterator);

  if (s.substr(0, iterator) == name_in_topo)
    return true;
  else
    return false;
}

bool utils::AtomSpecifier::_checkResName(int m, int a, std::string s)const
{
  if (m<0) throw Exception("checking residue name for solvent not possible");

  std::string::size_type iterator=s.find('?');
  std::string name_in_topo;

  const int resn = d_sys->mol(m).topology().resNum(a);
  name_in_topo=d_sys->mol(m).topology().resName(resn).substr(0, iterator);

  if (s.substr(0, iterator) == name_in_topo)
    return true;
  else
    return false;
}

int utils::AtomSpecifier::_expandSolvent()const
{
  int nsa = d_sys->sol(0).topology().numAtoms();
  d_nsm = d_sys->sol(0).numPos() / nsa;

  // first remove all atoms that are in the list due to an earlier
  // expansion. These have d_mol[i] == -2
  for(unsigned int i=0; i<d_specatom.size(); ++i){
    if(d_specatom[i]->mol()==-2) {
      _unexpandSolvent(i);
      i--;
    }
  }

  // now add the atoms for every molecule
  for(int i=0; i < d_nsm; ++i){
    for(unsigned int j=0; j< d_solventType.size(); ++j){
      _appendSolvent(-2, i*nsa+d_solventType[j]);
    }
  }

  return d_specatom.size();
}

bool utils::AtomSpecifier::_expand()const
{
  return (d_nsm !=
	  d_sys->sol(0).numPos() / d_sys->sol(0).topology().numAtoms());
}

utils::AtomSpecifier::AtomSpecifier(gcore::System &sys)
{
  d_sys = &sys;
  d_not_atoms = NULL;
  // set d_nsm to something weird so that it will be expanded on first use
  d_nsm = -1;
}

utils::AtomSpecifier::AtomSpecifier(gcore::System &sys, string s)
{
  d_sys = &sys;
  d_not_atoms = NULL;
  // set d_nsm to something weird so that it will be expanded on first use
  d_nsm = -1;
  parse(s);
}

utils::AtomSpecifier::AtomSpecifier(gcore::System &sys, string s, int x)
{
  d_sys = &sys;
  d_not_atoms = NULL;
  // set d_nsm to something weird so that it will be expanded on first use
  d_nsm = -1;
  parse(s, x);
}

utils::AtomSpecifier::~AtomSpecifier()
{
  // delete all atoms from the spec
  // (this is slightly undefined behaviour, as the deleted pointers
  //  are left in the vector (for clear))
  for(std::vector<SpecAtom *>::iterator it = d_specatom.begin(),
	to = d_specatom.end(); it != to; ++it){
    delete *it;
  }
  if (d_not_atoms != NULL)
    delete d_not_atoms;
}

void utils::AtomSpecifier::setSystem(gcore::System &sys)
{
  d_sys = &sys;
  // set the number of solvent molecules to something weird so that
  // we are sure the will be correctly expanded for the new system
  d_nsm = -1;

  if (d_not_atoms != NULL)
    d_not_atoms->setSystem(sys);

  // need to change the systems in the SpecAtom vector
  for(unsigned int i=0; i<d_specatom.size(); ++i)
    d_specatom[i]->setSystem(sys);
}

int utils::AtomSpecifier::addSpecifier(string s, int x)
{
  parse(s, x);
  return d_specatom.size();
}


int utils::AtomSpecifier::addSpecifierStrict(string s, int x)
{
  parseStrict(s, x);
  return d_specatom.size();
}

int utils::AtomSpecifier::appendSpecifier(string s, int x)
{
    AtomSpecifier atom(*d_sys);
    atom.addSpecifier(s,x);
    int m=atom.mol(0);
    int a=atom.atom(0);

    if(m >= d_sys->numMolecules())
        throw utils::AtomSpecifier::Exception(" molecule number out of range.\n");
    if(m >= 0){
        if(a >= d_sys->mol(m).topology().numAtoms())
            throw utils::AtomSpecifier::Exception(" atom number out of range.\n");
        d_specatom.push_back(new SpecAtom(*d_sys, m, a));
   }
    else{
        d_specatom.push_back(new SolventSpecAtom(*d_sys, m, a));
   }
    return d_specatom.size();
}

int utils::AtomSpecifier::addAtom(int m, int a) {
  if (d_sys != NULL) { // only check if we have a system!
    if (m >= d_sys->numMolecules())
      throw utils::AtomSpecifier::Exception(" molecule number out of range.\n");
    if (m >= 0)
      if (a >= d_sys->mol(m).topology().numAtoms())
        throw utils::AtomSpecifier::Exception(" atom number out of range.\n");
  }
  _appendAtom(m,a);

  return d_specatom.size();
}

int utils::AtomSpecifier::addGromosAtom(int a)
{
  int m=0;

  while(a >= d_sys->mol(m).numAtoms()){
    a -= d_sys->mol(m).numAtoms();
    m++;
    if(m >= d_sys->numMolecules()){
      m=-1;
      break;
    }
  }

  if(m >= 0)
    if(a >= d_sys->mol(m).topology().numAtoms())
      throw utils::AtomSpecifier::Exception(" atom number out of range.\n");

  _appendAtom(m,a);
  return d_specatom.size();
}

int utils::AtomSpecifier::addMolecule(int m)
{
  if(m >= int(d_sys->numMolecules()))
    throw utils::AtomSpecifier::Exception(" molecule number out of range.\n");

  for(int i=0; i < d_sys->mol(m).numAtoms(); ++i)
    _appendAtom(m, i);

  return d_specatom.size();

}

int utils::AtomSpecifier::addAtomStrict(int m, int a)
{
  if(m >= int(d_sys->numMolecules()))
    throw utils::AtomSpecifier::Exception(" molecule number out of range.\n");
  if(m >= 0){
    if(a >= d_sys->mol(m).topology().numAtoms())
      throw utils::AtomSpecifier::Exception(" atom number out of range.\n");

    d_specatom.push_back(new SpecAtom(*d_sys, m, a));

  }
  else{
    d_specatom.push_back(new SolventSpecAtom(*d_sys, m, a));
  }

  return d_specatom.size();
}

int utils::AtomSpecifier::addType(std::string s)
{
  //loop over all solute atoms
  for(int m = 0; m < d_sys->numMolecules(); ++m)
    addType(m, s);
  // and do solvent
  addType(-1, s);

  return d_specatom.size();
}


int utils::AtomSpecifier::addType(int m, std::string s)
{
  //loop over all atoms
  if(m<0)
    addSolventType(s);
  else{
    for(int j=0; j<d_sys->mol(m).numAtoms(); ++j)
      if(_checkName(m, j, s))
	_appendAtom(m,j);
  }
  return d_specatom.size();
}

int utils::AtomSpecifier::addType(int m, std::string s, int beg, int end)
{
  //loop over all atoms
  if(m<0)
    throw Exception("Solvent: addType for a range of atoms not possible");

  if (end > d_sys->mol(m).numAtoms())
    throw Exception("addType: end of range out of range");
  if (beg < 0)
    throw Exception("addType: begin of range < 0");

  for(int j=beg; j<end; ++j)
    if(_checkName(m, j, s))
      _appendAtom(m,j);

  return d_specatom.size();
}

int utils::AtomSpecifier::addTypeStrict(int m, std::string s)
{
  //loop over all atoms
  if(m<0)
    addSolventType(s);
  else{
    for(int j=0; j<d_sys->mol(m).numAtoms(); ++j)
      if(_checkName(m, j, s))
          addAtomStrict(m,j);
//	_appendAtom(m,j);
  }
  return d_specatom.size();
}

int utils::AtomSpecifier::addTypeStrict(int m, std::string s, int beg, int end)
{
  //loop over all atoms
  if(m<0)
    throw Exception("Solvent: addType for a range of atoms not possible");

  if (end > d_sys->mol(m).numAtoms())
    throw Exception("addType: end of range out of range");
  if (beg < 0)
    throw Exception("addType: begin of range < 0");

  for(int j=beg; j<end; ++j)
    if(_checkName(m, j, s))
        addAtomStrict(m,j);
//      _appendAtom(m,j);

  return d_specatom.size();
}

int utils::AtomSpecifier::addSolventType(std::string s)
{
  for(int j=0; j < d_sys->sol(0).topology().numAtoms(); ++j){
    if(_checkName(-1,j,s)){
      int found=0;
      for(unsigned int i = 0; i < d_solventType.size(); ++i)
	    if(d_solventType[i]==j) found = 1;
      if(!found)
	    d_solventType.push_back(j);
    }
  }
  _expandSolvent();
  return d_solventType.size();
}

int utils::AtomSpecifier::addSolventType() const //since there is only 1 solvent allowed, we dont need to specify it...
{
  for(int j=0; j < d_sys->sol(0).topology().numAtoms(); ++j){ //get the solvent from topology
    int found=0;
    for(unsigned int i = 0; i < d_solventType.size(); ++i)
      if(d_solventType[i]==j) found = 1;
    if(!found)
      d_solventType.push_back(j); //d_solventtype is mutable...
  }
//cout << "add solvtype: " << d_solventType.size() << endl;
  //_expandSolvent();
  return d_solventType.size();
}

bool utils::AtomSpecifier::_compare(int i, int m, int a)const
{
  if(d_specatom[i]->mol() == m) return d_specatom[i]->atom() > a;
  if(d_specatom[i]->mol() >= 0 && m >= 0) return d_specatom[i]->mol() > m;
  if(d_specatom[i]->mol() <  0 && m <  0) return d_specatom[i]->atom() > a;
  if(d_specatom[i]->mol() <  0 && m >= 0) return true;
  if(d_specatom[i]->mol() >= 0 && m <  0) return false;
  return false;
}

void utils::AtomSpecifier::sort()
{
  for(unsigned int i=1; i < d_specatom.size(); ++i){
    SpecAtom * t = d_specatom[i];
    int j=i-1;
    while ((j >= 0) && _compare(j, t->mol(), t->atom())){
      d_specatom[j+1] = d_specatom[j];
      j--;
    }
    d_specatom[j+1] = t;
  }
}

int utils::AtomSpecifier::removeAtom(int m, int a)
{
  int i=findAtom(m,a);
  return removeAtom(i);
}

int utils::AtomSpecifier::removeAtom(int i)
{
  if(i < int(d_specatom.size()) && i >= 0){
    vector<SpecAtom *>::iterator it=d_specatom.begin() + i;
    SpecAtom * d = *it;

    d_specatom.erase(it);
    delete d;
  }
  return d_specatom.size();
}

int utils::AtomSpecifier::_unexpandSolvent(int i)const
{
  if(i < int(d_specatom.size()) && i >= 0){
    vector<SpecAtom *>::iterator it=d_specatom.begin() + i;
    SpecAtom * d = *it;

    d_specatom.erase(it);
    delete d;
  }
  return d_specatom.size();
}

int utils::AtomSpecifier::findAtom(int m, int a)const
{
  vector<SpecAtom *>::iterator it=d_specatom.begin(),
    it_to = d_specatom.end();

  int counter=0;
  // a bit of a nuisance that m=-1 and m=-2 could both mean the same in this
  // function. So split up the cases
  if(m<0){
    for( ; it != it_to; ++it, ++counter)
      if((*it)->mol() < 0 && (*it)->atom() == a)
	return counter;
  }

  for( ; it != it_to; ++it, ++counter)
    if((*it)->mol() == m && (*it)->atom() == a)
      return counter;

  return -1;
}

/**
 * copy constructor.
 */
utils::AtomSpecifier::AtomSpecifier(utils::AtomSpecifier const & as)
{
  if (this != &as){
    d_not_atoms = NULL;
    clear();

    d_sys=as.d_sys;
    d_nsm=as.d_nsm;

    std::vector<SpecAtom *>::const_iterator
      it = as.d_specatom.begin(),
      to = as.d_specatom.end();
    for( ; it != to; ++it){
      d_specatom.push_back((*it)->clone());
    }

    if (as.d_not_atoms != NULL) {
      d_not_atoms = new AtomSpecifier(*(as.d_not_atoms));
    }
    d_solventType = as.d_solventType;
  }
}

utils::AtomSpecifier & utils::AtomSpecifier::operator=(const utils::AtomSpecifier &as)
{
  if(this != &as){

    clear();

    d_sys=as.d_sys;
    d_nsm=as.d_nsm;

    for(unsigned int i=0; i < as.d_specatom.size(); ++i){
      d_specatom.push_back(as.atom()[i]->clone());
    }
    if (d_not_atoms != NULL) {
      delete d_not_atoms;
    }
    if (as.d_not_atoms != NULL)
      d_not_atoms = new AtomSpecifier(*(as.d_not_atoms));
    else
      d_not_atoms = NULL;

    d_solventType = as.d_solventType;
  }

  return *this;
}

utils::AtomSpecifier utils::AtomSpecifier::operator+(const AtomSpecifier &as)
{
  utils::AtomSpecifier temp(*as.d_sys);
  temp = *this;

  for(unsigned int i = 0; i < as.d_specatom.size(); ++i)
    temp.addAtom(as.d_specatom[i]->mol(), as.d_specatom[i]->atom());

  for(unsigned int i=0;i<as.d_solventType.size(); i++)
    temp.addSolventType(d_sys->sol(0).topology().atom(as.d_solventType[i]).name());

  if (as.d_not_atoms != NULL) {
    if (d_not_atoms == NULL)
      d_not_atoms = new AtomSpecifier(*(as.d_not_atoms));
    else
      *d_not_atoms = *d_not_atoms + *(as.d_not_atoms);
  }

  // if copy construction works here, why not use it in the first statement???
  return temp;
}

gmath::Vec *utils::AtomSpecifier::coord(int i)
{
  if(_expand()) _expandSolvent();
  return &d_specatom[i]->pos();
}

gmath::Vec & utils::AtomSpecifier::pos(int i)
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->pos();
}
gmath::Vec const & utils::AtomSpecifier::pos(int i)const
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->pos();
}
gmath::Vec & utils::AtomSpecifier::cosDisplacement(int i)
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->cosDisplacement();
}
gmath::Vec const & utils::AtomSpecifier::cosDisplacement(int i)const
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->cosDisplacement();
}
gmath::Vec & utils::AtomSpecifier::vel(int i)
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->vel();
}
gmath::Vec const & utils::AtomSpecifier::vel(int i)const
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->vel();
}

std::string utils::AtomSpecifier::name(int i)const
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->name();
}

int utils::AtomSpecifier::iac(int i)const
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->iac();
}

double utils::AtomSpecifier::radius(int i)const
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->radius();
}

double utils::AtomSpecifier::charge(int i)const
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->charge();
}
bool utils::AtomSpecifier::isPolarisable(int i)const
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->isPolarisable();
}
double utils::AtomSpecifier::poloffsiteGamma(int i)const
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->poloffsiteGamma();
}
int utils::AtomSpecifier::poloffsiteI(int i)const
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->poloffsiteI();
}
int utils::AtomSpecifier::poloffsiteJ(int i)const
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->poloffsiteJ();
}
double utils::AtomSpecifier::cosCharge(int i)const
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->cosCharge();
}
double utils::AtomSpecifier::mass(int i)const
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->mass();
}

std::string utils::AtomSpecifier::resname(int i)const
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->resname();
}

int utils::AtomSpecifier::resnum(int i)const
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->resnum();
}

int utils::AtomSpecifier::mol(int i)const
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->mol();
}

int utils::AtomSpecifier::atom(int i)const
{
  if(_expand()) _expandSolvent();
  return d_specatom[i]->atom();
}

int utils::AtomSpecifier::gromosAtom(int i)const
{
  if(_expand()) _expandSolvent();
  int maxmol=d_specatom[i]->mol();
  if(maxmol<0) maxmol=d_sys->numMolecules();
  int grom=d_specatom[i]->atom();
  for(int j=0; j< maxmol; ++j) grom+=d_sys->mol(j).numAtoms();
  return grom;
}

bool utils::AtomSpecifier::empty()const
{
  return !(d_specatom.size() > 0 ||
	  d_solventType.size() > 0);
}

unsigned int utils::AtomSpecifier::size()const
{
  if(_expand()) _expandSolvent();
  return d_specatom.size();
}

void utils::AtomSpecifier::clear()
{
  std::vector<SpecAtom *>::iterator
    it = d_specatom.begin(),
    to = d_specatom.end();

  for( ; it != to; ++it)
    delete *it;

  d_specatom.clear();
  d_solventType.resize(0);
  d_nsm=-1;

  if (d_not_atoms != NULL)
    d_not_atoms->clear();
}

std::vector<std::string> utils::AtomSpecifier::toString()const
{
  std::vector<std::string> s;
  ostringstream os;
  int m = -999, a_last = -999;
  bool in_range = false, first = true;

  if(_expand()) _expandSolvent();

  // std::cerr << "toTitle" << std::endl;
  // std::cerr << "\tspecatoms " << d_specatom.size() << std::endl;

  for(unsigned int i = 0; i < d_specatom.size(); ++i){

    // virtual atom
    if (d_specatom[i]->type() == spec_virtual){
      if (in_range){
	in_range = false;
	os << a_last + 1;
      }
      os << "," << d_specatom[i]->toString();
      m = -999;
      a_last = -999;
      continue;
    }

    // new molecule?
    if (d_specatom[i]->mol() != m){
      // std::cerr << "\t" << i << " new molecule" << std::endl;
      m = d_specatom[i]->mol();

      if (in_range){
	in_range = false;
	os << "-" << a_last + 1;
      }

      a_last = d_specatom[i]->atom();

      if (!first){
	os << ";";
	if (m < 0) os << "s";
	else os << m + 1;

	os << ":" << a_last + 1;
      }
      else{
	if (m < 0) os << "s";
	else os << m+1;

	os << ":" << a_last + 1;
	first = false;
      }
    }
    else if (a_last == d_specatom[i]->atom()-1){
      // std::cerr << "\t" << i << " not new molecule but in order" << std::endl;
      a_last = d_specatom[i]->atom();
      in_range = true;
    }
    else if (in_range){
      // std::cerr << "\twere in range, but not anymore" << std::endl;
      in_range = false;
      os << "-" << a_last + 1 << ",";
      a_last = d_specatom[i]->atom();
      os << a_last + 1;
    }
    else{
      // std::cerr << "\tjust an other unconnected atom" << std::endl;
      a_last = d_specatom[i]->atom();
      os << "," << a_last + 1;
    }
  }

  if (in_range){
    os << "-" << a_last+1;
  }

  // but the why???
  s.push_back(os.str());
  return s;

}

std::string utils::AtomSpecifier::toString(int i)const
{
  if(_expand()) _expandSolvent();
  ostringstream os;
  // virtual atom
  if (d_specatom[i]->type() == spec_virtual){
      os  << d_specatom[i]->toString();
      return os.str();
  }
  else if(d_specatom[i]->mol() < 0) os << "s";
  else os << d_specatom[i]->mol()+1;
  os << ":" << d_specatom[i]->atom()+1;
  return os.str();
}

