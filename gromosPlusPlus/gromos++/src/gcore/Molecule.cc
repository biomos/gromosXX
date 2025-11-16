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

// gcore_Molecule.cc
/**
 * Class Molecule
 * Addition: velocity configuration added to Molecule definition;
 * Author  : gee          
 */

#include "Molecule.h"

#include <cassert>
#include <vector>

#include "MoleculeTopology.h"
#include "../gmath/Vec.h"

using namespace std;
using gcore::Molecule;
using gcore::MoleculeTopology;
using gmath::Vec;

Molecule::Molecule(const MoleculeTopology &mt):
  d_mt(new MoleculeTopology(mt)),
  d_pos(0),
  d_vel(0),
  d_bfac(0),
  d_cosDisplacement(0){}

Molecule::Molecule(const Molecule &mol):
  d_mt(new MoleculeTopology(*mol.d_mt)),
  d_pos(mol.d_pos.size()),
  d_vel(mol.d_vel.size()),
  d_bfac(mol.d_bfac.size()),
  d_cosDisplacement(mol.d_cosDisplacement.size())
{
  for(int i=0; i<mol.numPos();++i)
  {
    d_pos[i]=new Vec(mol.pos(i));
  }
  for(int i=0; i<mol.numVel();++i){
    d_vel[i]=new Vec(mol.vel(i));
  }
  for(int i=0; i<mol.numBfac();++i){
    d_bfac[i]=0.0;
  }
  for(int i=0; i<mol.numCosDisplacements();++i){
    d_cosDisplacement[i] = new Vec(mol.cosDisplacement(i));
  }
}

Molecule::~Molecule(){
  delete d_mt;
  for(unsigned int i=0; i<d_pos.size();++i)
    delete d_pos[i];

  for(unsigned int i=0; i<d_vel.size();++i)
    delete d_vel[i];

  for(unsigned int i=0; i<d_cosDisplacement.size();++i)
    delete d_cosDisplacement[i];
}

void Molecule::initPos(){
  d_pos.resize(numAtoms());
  for(int i=0; i < numAtoms(); ++i){
    d_pos[i]=new Vec();
  }
}
void Molecule::initVel(){
  d_vel.resize(numAtoms());
  for(int i=0; i < numAtoms(); ++i){
    d_vel[i]=new Vec();
  }
}
void Molecule::initBfac(){
  d_bfac.resize(numAtoms());
  for(int i=0; i < numAtoms(); ++i){
    d_bfac[i]=0.0;
  }
}
void Molecule::setBfac(int i, double b){
    d_bfac[i]=b;
}
void Molecule::initCosDisplacements(){
  d_cosDisplacement.resize(numAtoms());
  for(int i=0; i < numAtoms(); ++i){
    d_cosDisplacement[i]=new Vec();
  }
}
MoleculeTopology &Molecule::topology()
{
  return *d_mt;
}

const MoleculeTopology &Molecule::topology()const{
  return *d_mt;
}

int Molecule::numAtoms()const{
  return d_mt->numAtoms();
}

