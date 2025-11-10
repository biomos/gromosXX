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

// gcore_VirtualAtoms.cc
/**
 * Class VirtualAtoms
 */

#include <cassert>
#include <vector>

#include "System.h"
#include "VirtualAtoms.h"
#include "../gcore/Exclusion.h"
#include "../utils/VirtualAtom.h"
#include "../utils/AtomSpecifier.h"
#include "../gcore/GromosForceField.h"
#include "../gcore/VirtualAtomType.h"
#include "../gromos/Exception.h"

using namespace std;
using gcore::VirtualAtoms;

VirtualAtoms::VirtualAtoms(){}

VirtualAtoms::VirtualAtoms(utils::AtomSpecifier as, gcore::GromosForceField &gff)
{
  for(unsigned int i =0; i< as.size(); i++){
    if(as.atom()[i]->type()==utils::spec_virtual){ 
      utils::VirtualAtom va(*(as.sys()), 
 	                    as.atom()[i]->conf(), 
                            static_cast<utils::VirtualAtom::virtual_type>(as.atom()[i]->virtualType()));
      va.setDish(gff.virtualAtomType(as.atom()[i]->virtualType()).dis1());
      va.setDisc(gff.virtualAtomType(as.atom()[i]->virtualType()).dis2());
      d_vas.push_back(va);
      d_iac.push_back(-1);
      d_charge.push_back(0.0);
      d_exclusion.push_back(gcore::Exclusion());
      d_exclusion14.push_back(gcore::Exclusion());
    } else if(as.atom()[i]->type()==utils::spec_solute){
      utils::VirtualAtom va(*(as.sys()),
                            as.mol(i),
                            as.atom(i),
                            static_cast<utils::VirtualAtom::virtual_type>(0));
      va.setDish(gff.virtualAtomType(0).dis1());
      va.setDisc(gff.virtualAtomType(0).dis2());
      d_vas.push_back(va);
      d_iac.push_back(as.iac(i));
      d_charge.push_back(as.charge(i));
      d_exclusion.push_back(as.sys()->mol(as.mol(i)).topology().atom(as.atom(i)).exclusion());
      d_exclusion14.push_back(as.sys()->mol(as.mol(i)).topology().atom(as.atom(i)).exclusion14());
    } else {
      throw gromos::Exception("VirtualAtoms", "Don't know how to handle atomspecifier");
    }
  }
}
    
VirtualAtoms::VirtualAtoms(const VirtualAtoms &vas):
   d_vas(vas.d_vas),
   d_iac(vas.d_iac),
   d_charge(vas.d_charge),
   d_exclusion(vas.d_exclusion),
   d_exclusion14(vas.d_exclusion14)
{
}

VirtualAtoms::~VirtualAtoms()
{
  d_vas.erase(d_vas.begin(), d_vas.end());
}
void VirtualAtoms::addVirtualAtom(utils::AtomSpecifier as, gcore::GromosForceField &gff, int iac, double charge, gcore::Exclusion e, gcore::Exclusion e14)
{
for(unsigned int i =0; i< as.size(); i++){
    if(as.atom()[i]->type()==utils::spec_virtual){
      utils::VirtualAtom va(*(as.sys()),
                            as.atom()[i]->conf(),
                            static_cast<utils::VirtualAtom::virtual_type>(as.atom()[i]->virtualType()));
      va.setDish(gff.virtualAtomType(as.atom()[i]->virtualType()).dis1());
      va.setDisc(gff.virtualAtomType(as.atom()[i]->virtualType()).dis2());
      d_vas.push_back(va);
      d_iac.push_back(iac);
      d_charge.push_back(charge);
      d_exclusion.push_back(e);
      d_exclusion14.push_back(e14);
    } else if(as.atom()[i]->type()==utils::spec_solute){
      utils::VirtualAtom va(*(as.sys()),
                            as.mol(i),
                            as.atom(i),
                            static_cast<utils::VirtualAtom::virtual_type>(0));
      va.setDish(gff.virtualAtomType(0).dis1());
      va.setDisc(gff.virtualAtomType(0).dis2());
      d_vas.push_back(va);
      d_iac.push_back(as.iac(i));
      d_charge.push_back(as.charge(i));
      d_exclusion.push_back(as.sys()->mol(as.mol(i)).topology().atom(as.atom(i)).exclusion());
      d_exclusion14.push_back(as.sys()->mol(as.mol(i)).topology().atom(as.atom(i)).exclusion14());
    } else {
      throw gromos::Exception("VirtualAtoms", "Don't know how to handle atomspecifier");
    }
  }
}


void VirtualAtoms::addVirtualAtom(System &sys, std::vector<int> conf, int type,
                                  double dish, double disc, 
                                  int iac, double charge,
                                  gcore::Exclusion e, gcore::Exclusion e14)
{
  utils::VirtualAtom va(sys, static_cast<utils::VirtualAtom::virtual_type>(type), conf, dish, disc, 0);
  d_vas.push_back(va);
  d_iac.push_back(iac);
  d_charge.push_back(charge);
  d_exclusion.push_back(e);
  d_exclusion14.push_back(e14);
}
unsigned int VirtualAtoms::numVirtualAtoms()const{
     return d_vas.size();
}
utils::VirtualAtom VirtualAtoms::atom(int i){
    assert(i < d_vas.size());
    return d_vas[i];
}
const utils::VirtualAtom VirtualAtoms::atom(int i)const{
    assert(i < d_vas.size());
    return d_vas[i];
}
void VirtualAtoms::setSystem(gcore::System &sys){
  for(unsigned int i=0; i< d_vas.size(); i++){
    d_vas[i].setSystem(sys);
  }
}
