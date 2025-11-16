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

// gcore_BbLink.cc

#include "BbLink.h"

#include <cassert>

#include "MoleculeTopology.h"

using namespace std;
using namespace gcore;

BbLink::BbLink(const BbLink& mt)
{
  for(int i=0; i<mt.numAtoms(); i++)
    MoleculeTopology::addAtom(mt.atom(i));
  
  for(int i=0; i<mt.numPexcl(); i++)
    addPexcl(mt.pexcl(i));
  
  BondIterator bi(mt);
  for(;bi;++bi)
    MoleculeTopology::addBond(bi());
  
  BondDipoleIterator bdi(mt);
  for(;bdi;++bdi)
    MoleculeTopology::addDipoleBond(bdi());
  
  AngleIterator ai(mt);
  for(;ai;++ai)
    MoleculeTopology::addAngle(ai());
  
  DihedralIterator di(mt);
  for(;di;++di)
    MoleculeTopology::addDihedral(di());
  
  ImproperIterator ii(mt);
  for(;ii;++ii)
    MoleculeTopology::addImproper(ii());
  
  LJExceptionIterator lji(mt);
  for(;lji;++lji)
    MoleculeTopology::addLJException(lji());
  
  setResName(mt.resName());
  setRep(mt.rep());

  for(int i=0; i < mt.numAtoms(); i++)
    setLinkRes(i, mt.linkRes(i));
}
// Methods

BbLink &BbLink::operator=(const BbLink &mt){
  if (this != &mt){
    this->BbLink::~BbLink();
    new(this) BbLink(mt);
  }
  return *this;
}

void BbLink::setLinkRes(const unsigned int a, const unsigned int i){
  if(d_linkres.size() <= a)
   d_linkres.resize(a+1); 
  d_linkres[a]=i;
}
  
int BbLink::linkRes(const int a)const{
  return d_linkres[a];
}
