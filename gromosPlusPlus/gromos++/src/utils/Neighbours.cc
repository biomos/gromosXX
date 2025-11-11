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
#include "Neighbours.h"

#include <cassert>
#include <vector>

#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/Solvent.h"
#include "../gcore/Bond.h"
#include "../gcore/Constraint.h"

using namespace gcore;
using namespace std;
using utils::Neighbours;

static void init(Neighbours *n, const MoleculeTopology &mt, int i){
  BondIterator iter(mt);
  for(;iter;++iter){
    if(iter()[0]==i)n->push_back(iter()[1]);
    if(iter()[1]==i)n->push_back(iter()[0]);
  }
}

static void inits(Neighbours *n, const SolventTopology &mt, int i){
  ConstraintIterator iter(mt);
  for(;iter;++iter){
    if(iter()[0]==i)n->push_back(iter()[1]);
    if(iter()[1]==i)n->push_back(iter()[0]);
  }
}

Neighbours::Neighbours(const System &sys, int mol, int i):
  vector<int>()
{
  //  int j=0, atNum=0;
  //  while(i > (atNum+=sys.mol(j).numAtoms())) ++j;
  //  int k=i-atNum+sys.mol(j).numAtoms();

  init(this,sys.mol(mol).topology(),i);
}

Neighbours::Neighbours(const System &sys, int mol, int i, int j):
  vector<int>()
{
  //  int j=0, atNum=0;
  //  while(i > (atNum+=sys.mol(j).numAtoms())) ++j;
  //  int k=i-atNum+sys.mol(j).numAtoms();

  inits(this,sys.sol(mol).topology(),i);
}

Neighbours::Neighbours(const MoleculeTopology &mt, int i): vector<int>()
{
  init(this,mt,i);
}

Neighbours::Neighbours(const Molecule &mol, int k):
  vector<int>()
{
  init(this,mol.topology(),k);
}





