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

// utils_Rmsd.cc
#include "Rmsd.h"

#include <cassert>
#include <vector>

#include "../fit/Reference.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gmath/Vec.h"
#include "../utils/AtomSpecifier.h"
#include "../utils/Value.h"
#include "../utils/Property.h"
#include "../utils/PropertyContainer.h"

using gcore::System;
using utils::Rmsd;
using fit::Reference;

using namespace std;

Rmsd::Rmsd(const Reference *ref){
  d_ref=ref;
}



double Rmsd::rmsd(const System &sys){
  double rmsd2=0;
  for(int m=0;m<sys.numMolecules();++m)
    for(int n=0;n<sys.mol(m).numAtoms();++n)
      if(d_ref->weight(m,n))
	rmsd2+=d_ref->weight(m,n)*
	  (d_ref->sys().mol(m).pos(n) 
	   - sys.mol(m).pos(n)).abs2();

  return sqrt(rmsd2);
}


/**
 * @todo should include periodicity. But this depends on the property?
 */
double Rmsd::rmsdproperty(const System &sys){

  double rmsd2=0;
  PropertyContainer prop_ref = *d_prop_ref;
  PropertyContainer prop_sys = *d_prop_sys;

  for(unsigned int i=0; i < d_prop_ref->size(); i++){
    
    Value res = abs2(prop_ref[i]->calc() - prop_sys[i]->calc());
    rmsd2 += res.scalar();
  }
  
  return sqrt(rmsd2/d_prop_ref->size());
}

void Rmsd::addproperty(const PropertyContainer *prop_ref, const PropertyContainer *prop_sys) {  
  d_prop_ref = prop_ref;
  d_prop_sys = prop_sys;
}
