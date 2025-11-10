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

// GromosForceField.cc

#include "GromosForceField.h"

#include <vector>
#include <map>
#include <string>

#include "AtomPair.h"
#include "MassType.h"
#include "BondType.h"
#include "AngleType.h"
#include "DihedralType.h"
#include "ImproperType.h"
#include "LJExceptionType.h"
#include "LJType.h"
#include "CGType.h"
#include "VirtualAtomType.h"
#include "../gromos/Exception.h"


using namespace std;

namespace gcore{

class GromosForceField_i{
  friend class GromosForceField;
  double d_fpepsi, d_hbar, d_spdl, d_boltz;
  std::string d_ffcode;
  vector<string> d_atomTypeName;
  map<int, MassType> d_massType;
  map<int, VirtualAtomType> d_virtualAtomType;
  map<int, BondType> d_bondType;
  map<int, AngleType> d_angleType;
  map<int, DihedralType> d_dihedralType;
  map<int, ImproperType> d_improperType;
  // this map is needed to remember the LJ exception types read from the parameter file
  map<int, LJExceptionType> d_ljExceptionType;
  // and this list is filled additionally, in case a topology is read
  map<AtomPair,LJExceptionType> d_ljException;
  map<AtomPair,LJType> d_ljType;
  map<AtomPair,CGType> d_cgType;
  GromosForceField_i():
    d_fpepsi(0), d_hbar(0), d_spdl(0), d_boltz(0), d_ffcode("_no_FORCEFIELD_block_given_"),
    d_atomTypeName(), d_massType(), d_virtualAtomType(), d_bondType(), d_angleType(),
    d_dihedralType(), d_improperType(), d_ljExceptionType(), d_ljType(), d_cgType()
  {}
  GromosForceField_i(const GromosForceField_i &gff):
    d_fpepsi(gff.d_fpepsi), d_hbar(gff.d_hbar), d_spdl(gff.d_spdl),
    d_boltz(gff.d_boltz), d_ffcode(gff.d_ffcode),
    d_atomTypeName(gff.d_atomTypeName), d_massType(gff.d_massType),
    d_virtualAtomType(gff.d_virtualAtomType),
    d_bondType(gff.d_bondType), d_angleType(gff.d_angleType),
    d_dihedralType(gff.d_dihedralType), d_improperType(gff.d_improperType),
    d_ljExceptionType(gff.d_ljExceptionType), d_ljType(gff.d_ljType),
    d_cgType(gff.d_cgType)
  {}

  ~GromosForceField_i(){}
};

GromosForceField::GromosForceField(): d_this(new GromosForceField_i()){}

  GromosForceField::GromosForceField ( const GromosForceField &gff):
    d_this(new GromosForceField_i(*gff.d_this)){}


GromosForceField::~GromosForceField(){delete d_this;}

void GromosForceField::setFpepsi(double fpepsi)
{d_this->d_fpepsi=fpepsi;}

void GromosForceField::setHbar(double hbar)
{d_this->d_hbar=hbar;}

void GromosForceField::setSpdl(double spdl)
{d_this->d_spdl=spdl;}

void GromosForceField::setBoltz(double boltz)
{d_this->d_boltz=boltz;}

void GromosForceField::setForceField(const string str)
{d_this->d_ffcode=str;}

void GromosForceField::addAtomTypeName(const string &str)
{d_this->d_atomTypeName.push_back(str);}

void GromosForceField::addMassType(const MassType &b)
{d_this->d_massType[b.n()] = b;}

void GromosForceField::addVirtualAtomType(const VirtualAtomType &va)
{d_this->d_virtualAtomType[va.code()] = va;}

void GromosForceField::addBondType(const BondType &b)
{d_this->d_bondType[b.code()] = b;}

void GromosForceField::addAngleType(const AngleType &b)
{d_this->d_angleType[b.code()] = b;}

void GromosForceField::addDihedralType(const DihedralType &b)
{d_this->d_dihedralType[b.code()] = b;}

void GromosForceField::addImproperType(const ImproperType &b)
{d_this->d_improperType[b.code()] = b;}

void GromosForceField::addLJExceptionType(const LJExceptionType &b)
{d_this->d_ljExceptionType[b.code()] = b;}

void GromosForceField::setLJException(const AtomPair &p, const LJExceptionType &l)
{ d_this->d_ljException[p]=l;}

void GromosForceField::setLJType(const AtomPair &p, const LJType &l)
{ d_this->d_ljType[p]=l;}

void GromosForceField::setCGType(const AtomPair &p, const CGType &l)
{ d_this->d_cgType[p]=l;}

double GromosForceField::fpepsi()const{
  return d_this->d_fpepsi;}

double GromosForceField::hbar()const{
  return d_this->d_hbar;}

double GromosForceField::spdl() const{
  return d_this->d_spdl;}

double GromosForceField::boltz()const{
  return d_this->d_boltz;}

std::string GromosForceField::ForceField()const{
  return d_this->d_ffcode;}

int GromosForceField::numAtomTypeNames()const
{ return d_this->d_atomTypeName.size();}

int GromosForceField::numMassTypes()const
{return d_this->d_massType.size();}

int GromosForceField::numVirtualAtomTypes()const
{ return d_this->d_virtualAtomType.size();}

int GromosForceField::numBondTypes()const
{ return d_this->d_bondType.size();}

int GromosForceField::numAngleTypes()const
{ return d_this->d_angleType.size();}

int GromosForceField::numImproperTypes()const
{ return d_this->d_improperType.size();}

int GromosForceField::numLJExceptionTypes()const
{ return d_this->d_ljExceptionType.size();}

int GromosForceField::numDihedralTypes()const
{ return d_this->d_dihedralType.size();}

int GromosForceField::numLJTypes()const
{ return d_this->d_ljType.size();}

int GromosForceField::numCGTypes()const
{ return d_this->d_cgType.size();}

const MassType &GromosForceField::massType(const int i) const
{return d_this->d_massType[i];}

double GromosForceField::findMass(const int i)const
{
  for(unsigned int k=0; k<d_this->d_massType.size(); k++)
    if(d_this->d_massType[k].n()==i) return d_this->d_massType[k].am();
  return 0.0;
}

int GromosForceField::findMassType(double mass) const
{
  for(unsigned int k = 0; k < d_this->d_massType.size(); ++k)
    if (d_this->d_massType[k].am() == mass) return int(d_this->d_massType[k].n());
  return -1;
}

const VirtualAtomType &GromosForceField::virtualAtomType(const int i) const
{ return d_this->d_virtualAtomType[i];}

const VirtualAtomType &GromosForceField::virtualAtomTypeLine(const int i) const
{ std::map<int,VirtualAtomType>::const_iterator iter = d_this->d_virtualAtomType.begin();
  for(int j=0; j< i; j++, ++iter){
    if (iter == d_this->d_virtualAtomType.end()){
      throw gromos::Exception("GromosForceField", "Trying to access a Virtual Atom Type that does not exist");
    }
  }
  return iter->second;
}

const bool GromosForceField::findVirtualAtomType(const int i) const
{
  for(int j=0; j< numVirtualAtomTypes(); j++)
    if (virtualAtomTypeLine(j).code() == i) return true;
  return false;
}
const BondType &GromosForceField::bondType(const int i) const
{ return d_this->d_bondType[i];}

const AngleType &GromosForceField::angleType(const int i) const
{ return d_this->d_angleType[i];}

const DihedralType &GromosForceField::dihedralType(const int i) const
{ return d_this->d_dihedralType[i];}

const ImproperType &GromosForceField::improperType(const int i) const
{ return d_this->d_improperType[i];}

const LJExceptionType &GromosForceField::ljExceptionType(const int i) const
{ return d_this->d_ljExceptionType[i];}

const map<AtomPair, LJExceptionType> &GromosForceField::ljException() const
{ return d_this->d_ljException;}

const string &GromosForceField::atomTypeName(const int i) const 
{ return d_this->d_atomTypeName[i];}

const LJType &GromosForceField::ljType(const AtomPair &p) const
{return d_this->d_ljType[p]; }

const CGType &GromosForceField::cgType(const AtomPair &p) const
{return d_this->d_cgType[p];
  }

  int GromosForceField::dummyAtomType() const {
    // first try to find the atom type named DUM
    int dummyAtomType = -1;
    for (int i = 0; i < numAtomTypeNames(); ++i) {
      if (atomTypeName(i) == "DUM")
        dummyAtomType = i;
    }

    if (dummyAtomType == -1)
      return -1;

    // then check whether this is really a dummy
    bool isDummy = true;
    for (int i = 0; i < numAngleTypes(); ++i) {
      const LJType lj = ljType(AtomPair(i, dummyAtomType));
      if (lj.c12() != 0.0 && lj.c6() != 0.0 && lj.cs12() != 0.0 && lj.cs6() != 0.0) {
        isDummy = false;
        break;
      } // if is dummy
    } // for atoms

    if (isDummy)
      return dummyAtomType;

    return -1;
  }

}
