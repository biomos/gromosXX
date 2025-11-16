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

// gcore_SolventTopology.cc

#include "SolventTopology.h"

#include <cassert>
#include <set>
#include <string>
#include <vector>

#include "AtomTopology.h"
#include "Constraint.h"

using namespace std;
using gcore::SolventTopology;
using gcore::ConstraintIterator;
using gcore::Constraint;
using gcore::AtomTopology;

class gcore::SolventTopology_i{

  friend class gcore::SolventTopology;
  friend class gcore::ConstraintIterator;

  vector<AtomTopology> d_atoms;
  set<Constraint> d_constraints;
  string d_solvName;
  
  SolventTopology_i():
    d_atoms(),
    d_constraints(),
    d_solvName()
  {d_solvName = "SOLV";}
  ~SolventTopology_i(){}
};


SolventTopology::SolventTopology() : 
  d_this(new SolventTopology_i())
{}

SolventTopology::SolventTopology(const SolventTopology& mt):
  d_this(new SolventTopology_i())
{
  d_this->d_atoms=(mt.d_this->d_atoms);
  d_this->d_constraints=(mt.d_this->d_constraints);
  d_this->d_solvName=(mt.d_this->d_solvName);
}

SolventTopology::~SolventTopology(){delete d_this;}

// Methods

SolventTopology &SolventTopology::operator=(const SolventTopology &mt){
  if (this != &mt){
    this->SolventTopology::~SolventTopology();
    new(this) SolventTopology(mt);
  }
  return *this;
}

void SolventTopology::addAtom(const AtomTopology &a){
  d_this->d_atoms.push_back(a);
  return;
}

void SolventTopology::addConstraint(const Constraint &b){
  // add checks if Constraint there?
  d_this->d_constraints.insert(b);
}

void SolventTopology::setSolvName(const string &s){

  d_this->d_solvName=s;
}

void SolventTopology::clearH()
{
  for(unsigned int i=0;i < d_this->d_atoms.size(); i++)
    d_this->d_atoms[i].setH(false);
}

void SolventTopology::setHmass(double mass)
{
  for(unsigned int i=0; i< d_this->d_atoms.size(); i++)
    if(d_this->d_atoms[i].mass() == mass)
      d_this->d_atoms[i].setH(true);
}

void SolventTopology::setHiac(int iac)
{
   for(unsigned int i=0; i< d_this->d_atoms.size(); i++)
    if(d_this->d_atoms[i].iac() == iac)
      d_this->d_atoms[i].setH(true);
}

int SolventTopology::numAtoms()const{return d_this->d_atoms.size();}

const AtomTopology &SolventTopology::atom(int i)const{
  assert(i < int(d_this->d_atoms.size()));
  return d_this->d_atoms[i];
}
AtomTopology &SolventTopology::atom(int i)
{
  assert(i < int(d_this->d_atoms.size()));
  return d_this->d_atoms[i];
}

const string &SolventTopology::solvName()const{
    return d_this->d_solvName;
}

int SolventTopology::numConstraints()const{return d_this->d_constraints.size();}

class gcore::ConstraintIterator_i{
  friend class gcore::ConstraintIterator;
  set<Constraint>::iterator d_it;
  const SolventTopology *d_mt;
  // not implemented
  ConstraintIterator_i(const ConstraintIterator_i&);
  ConstraintIterator_i &operator=(const ConstraintIterator_i &);
public:
  ConstraintIterator_i():
    d_it(){d_mt=0;}
};

ConstraintIterator::ConstraintIterator(const SolventTopology &mt):
  d_this(new ConstraintIterator_i())
{
  d_this->d_it=mt.d_this->d_constraints.begin();
  d_this->d_mt=&mt;
}

ConstraintIterator::~ConstraintIterator(){delete d_this;}

void ConstraintIterator::operator++(){
  ++(d_this->d_it);
}

const Constraint &ConstraintIterator::operator()()const{
  return *(d_this->d_it);
}

ConstraintIterator::operator bool()const{
  return d_this->d_it != d_this->d_mt->d_this->d_constraints.end();
}





