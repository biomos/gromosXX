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

// gcore_BbSolute.cc

#include "BbSolute.h"

#include <cassert>
#include <string>

#include "AtomTopology.h"
#include "Exclusion.h"
#include "MoleculeTopology.h"
#include "Exclusion.h"

using namespace std;
using namespace gcore;

BbSolute::BbSolute(const BbSolute& mt)
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
}

BbSolute::BbSolute(const MoleculeTopology & mt)
{
  for(int i=0; i<mt.numAtoms(); i++)
    MoleculeTopology::addAtom(mt.atom(i));
  
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
  
  setResName(mt.resName(0));
  setRep(0);
}

// Methods

BbSolute &BbSolute::operator=(const BbSolute &mt){
  if (this != &mt){
    this->BbSolute::~BbSolute();
    new(this) BbSolute(mt);
  }
  return *this;
}

void BbSolute::addPexcl(const Exclusion &a)
{
  d_pexcl.push_back(a);
  return;
}

void BbSolute::setResName(const string &s){
  MoleculeTopology::setResName(0,s);
}

int BbSolute::numPexcl()const{return d_pexcl.size();}

const Exclusion &BbSolute::pexcl(int i)const
{
  assert(i < int(d_pexcl.size()));
  return d_pexcl[i];
}

const string &BbSolute::resName()const{
  return MoleculeTopology::resName(0);
}
void BbSolute::setRep(int i)
{
  d_rep=i;
}

int BbSolute::rep()const{return d_rep;} 
