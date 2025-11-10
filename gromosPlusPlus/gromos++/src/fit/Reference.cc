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

// fit_Reference.cc

#include "Reference.h"

#include <cassert>
#include <vector>
#include <string>

#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/AtomTopology.h"
#include "../utils/AtomSpecifier.h"

using namespace gcore;
using namespace std;
using namespace utils;
using fit::Reference;

class fit::Reference_i{
  friend class fit::Reference;
  System *d_sys;
  vector<vector<double> > d_weights;
  Reference_i():d_weights(){}
  ~Reference_i(){}
};

Reference::Reference(gcore::System *sys):
  d_this(new Reference_i())
{
  d_this->d_sys=sys;
  d_this->d_weights.resize(sys->numMolecules());
  for(int i=0;i<sys->numMolecules();++i){
    d_this->d_weights[i].resize(sys->mol(i).numAtoms(),0);
  }
}

Reference::~Reference(){
  delete d_this;
}

void Reference::addClass(int mol, const string &name){
  for(int i=0;i<d_this->d_sys->mol(mol).numAtoms();++i){
    if(d_this->d_sys->mol(mol).topology().atom(i).name()==name || name=="ALL")
      d_this->d_weights[mol][i]=1;
  }
  rescale();
}

void Reference::addAtom(int m, int i){
  assert(m<int(d_this->d_weights.size()));
  assert(i<int(d_this->d_weights[m].size()));
  d_this->d_weights[m][i]=1;
  rescale();
}

void Reference::addAtomSpecifier(utils::AtomSpecifier as) {
  for (unsigned int i=0; i < as.size(); ++i) addAtom(as.mol(i), as.atom(i));
} 


void Reference::rescale(){
  // gives equal weight to all non 0 elements.
  int tot=0;
  const System &sys=*d_this->d_sys;
  for(int i=0;i<sys.numMolecules();++i)
    for(int j=0;j<sys.mol(i).numAtoms();++j)
      if(d_this->d_weights[i][j])tot++;
  double w=1.0/double(tot);
  for(int i=0;i<sys.numMolecules();++i)
    for(int j=0;j<sys.mol(i).numAtoms();++j)
      if(d_this->d_weights[i][j])d_this->d_weights[i][j]=w;
}

void Reference::setWeight(int m, int i, double w){
  assert(m<int(d_this->d_weights.size()));
  assert(i<int(d_this->d_weights[m].size()));
  d_this->d_weights[m][i]=w;
}

void Reference::normalise(){
  double tot=0;
  const System &sys=*d_this->d_sys;
  for(int i=0;i<sys.numMolecules();++i)
    for(int j=0;j<sys.mol(i).numAtoms();++j)
      tot+=d_this->d_weights[i][j];
  for(int i=0;i<sys.numMolecules();++i)
    for(int j=0;j<sys.mol(i).numAtoms();++j)
      d_this->d_weights[i][j]/=tot;
}  

 void Reference::makePosList(const System &sys,int molecule, const std::string &atomtype, std::vector<int> &poslist){
     if (atomtype=="ALL"){
       for (int j=0;j<sys.mol(molecule).numAtoms();++j)
       poslist.push_back(j);    
     }
     else {
      for (int j=0;j<sys.mol(molecule).numAtoms();++j)
       if (sys.mol(molecule).topology().atom(j).name()==atomtype){
       poslist.push_back(j);
       }
      }
}

const gcore::System &Reference::sys()const{
  return *d_this->d_sys;
}

gcore::System &Reference::sys(){
  return *d_this->d_sys;
}

double Reference::weight(int m, int i)const{
  return d_this->d_weights[m][i];
}
