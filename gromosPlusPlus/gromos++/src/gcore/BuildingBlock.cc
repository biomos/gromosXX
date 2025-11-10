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

// gcore_BuildingBlock.cc

#include "BuildingBlock.h"

#include <cassert>
#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <iostream>

#include "BbSolute.h"
#include "SolventTopology.h"
#include "../gromos/Exception.h"
#include "../args/Arguments.h"

using namespace gcore;

BuildingBlock::BuildingBlock():
  d_bb(),
  d_be(),
  d_bs(),
  d_fpepsi(),
  d_hbar(),
  d_spdl(),
  d_boltz(),
  d_physConstRead(false),
  //d_ffcode, is an empty set
  d_linkExclusions(),
  d_empty(true)
{}


BuildingBlock::BuildingBlock(const BuildingBlock &bld):
  d_bb(bld.d_bb.size()),
  d_be(bld.d_be.size()),
  d_bs(bld.d_bs.size()),
  d_fpepsi(bld.Fpepsi()),
  d_hbar(bld.Hbar()),
  d_spdl(bld.Spdl()),
  d_boltz(bld.Boltz()),
  d_physConstRead(bld.d_physConstRead),
  d_ffcode(bld.ForceField()),
  d_linkExclusions(bld.LinkExclusions()),
  d_empty(false)
{
  for (unsigned int i=0; i<d_bb.size();++i){
    d_bb[i]= new BbSolute(bld.bb(i));
  }
  for (unsigned int i=0; i<d_be.size();++i){
      d_be[i] = new BbSolute(bld.be(i));
  }
  for (unsigned int i=0; i<d_bs.size();++i){
      d_bs[i] = new SolventTopology(bld.bs(i));
  }
}

BuildingBlock::~BuildingBlock(){
  for (unsigned int i=0; i<d_bb.size();++i){
    delete d_bb[i];
  }
  for (unsigned int i=0; i<d_bs.size();++i){
    delete d_bs[i];
  }
  for (unsigned int i=0; i<d_be.size();++i){
    delete d_be[i];
  }
}

void BuildingBlock::addBuildingBlock(const BuildingBlock &bld)
{
  // check force field code, fpepsi, hbar and spdl
  if(!d_empty){
    if (bld.physConstRead()) {
      if (d_fpepsi != bld.Fpepsi()) {
        std::ostringstream os;
        os.precision(15);
        os.unsetf(std::ios::floatfield);
        os << "Value of FPEPSI is not identical in building block files.\n"
                << "Are there multiple non-identical PHYSICALCONSTANTS blocks in your mtb files?";
        throw gromos::Exception("BuildingBlock", os.str());
      }
      if (d_hbar != bld.Hbar()) {
        std::ostringstream os;
        os << "Value of HBAR is not identical in building block files.\n"
                << "Are there multiple non-identical PHYSICALCONSTANTS blocks in your mtb files?";
        throw gromos::Exception("BuildingBlock", os.str());
      }
      if (d_spdl != bld.Spdl()) {
        std::ostringstream os;
        os << "Value of SPDL is not identical in building block files.\n"
                << "Are there multiple non-identical PHYSICALCONSTANTS blocks in your mtb files?";
        throw gromos::Exception("BuildingBlock", os.str());
      }
      if (!(args::Arguments::inG96) && (d_boltz != bld.Boltz())) {
        std::ostringstream os;
        os << "Value of boltzmann constant is not identical in building block files.\n"
                << "Are there multiple non-identical PHYSICALCONSTANTS blocks in your mtb files?";
        throw gromos::Exception("BuildingBlock", os.str());
      }
    }
  } else{
    if(!bld.physConstRead()) {
      throw gromos::Exception("BuildingBlock", "No PHYSICALCONSTANTS block found "
              "in the main/first mtb file. Please make sure the first mtb file"
              " has a PHYSICALCONSTANTS block.");
    } else {
      d_fpepsi = bld.Fpepsi();
      d_hbar = bld.Hbar();
      d_spdl = bld.Spdl();
      d_boltz = bld.Boltz();
      d_ffcode = bld.ForceField();
      d_linkExclusions = bld.LinkExclusions();
      d_empty = false;
    }
  }
  
  // now add all the individual building blocks
  for(int i=0; i< bld.numBbSolutes(); ++i)
    addBbSolute(bld.bb(i));
  for(int i=0; i< bld.numBbEnds(); ++i)
    addBbEnd(bld.be(i));
  for(int i=0; i< bld.numBbSolvents(); ++i)
    addBbSolvent(bld.bs(i));
}

BuildingBlock &BuildingBlock::operator=(const BuildingBlock &bld){
  if(this != &bld){
    delete this;
    new(this) BuildingBlock(bld);
  }
  return *this;
}

void BuildingBlock::addBbSolute(const BbSolute &bb){
  d_bb.push_back(new BbSolute(bb));
}

void BuildingBlock::addBbSolvent(const SolventTopology &bs){
    d_bs.push_back(new SolventTopology(bs));
}

void BuildingBlock::addBbEnd(const BbSolute &be){
    d_be.push_back(new BbSolute(be));
}

int BuildingBlock::findBb(std::string s)
{
  //loop over all solute building blocks
  for(unsigned int i=0; i<d_bb.size(); ++i){
    if(d_bb[i]->resName()==s) return i+1;
  }
  //or, if not found over the end-building blocks
  for(unsigned int i=0; i<d_be.size(); ++i){
    if(d_be[i]->resName()==s) return -i-1;
  }
  return 0;
  
}
int BuildingBlock::findBb(std::string s, int &n)
{
  n=0;
  int index = 0;
  //loop over all solute building blocks
  for(unsigned int i=0; i<d_bb.size(); ++i){
    if(d_bb[i]->resName()==s) {n++; if(!index) index = i+1;}
    
  }
  //or, if not found over the end-building blocks
  for(unsigned int i=0; i<d_be.size(); ++i){
    if(d_be[i]->resName()==s) {n++; if(!index) index = -i-1; }
  }
  return index;
  
}
int BuildingBlock::findBs(std::string s)
{
  for(unsigned int i=0; i<d_bs.size(); ++i){
    if(d_bs[i]->solvName()==s) return i+1;
  }
  return 0;
}
int BuildingBlock::findBs(std::string s, int &n)
{
  n = 0;
  int index = 0;
  for(unsigned int i=0; i<d_bs.size(); ++i){
    if(d_bs[i]->solvName()==s) {
      n++;
      if(!index) {
        index = i + 1;
      }
    }
  }
  return index;
}







