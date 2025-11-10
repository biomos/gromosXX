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

// fit_TranslationalFit.cc

#include "TranslationalFit.h"

#include <cassert>

#include "Reference.h"
#include "PositionUtils.h"
#include "../gcore/System.h"
#include "../gmath/Vec.h"

using fit::TranslationalFit;
using fit::Reference;

// static Vec com(const System &sys, const Reference &ref);
// static Vec cog(const System &sys, const Reference &ref);


class fit::TranslationalFit_i{
  friend class fit::TranslationalFit;
  Reference *d_ref;
  TranslationalFit_i(Reference *ref)
  {
    d_ref=ref;
  }
  ~TranslationalFit_i(){}
};

TranslationalFit::TranslationalFit(Reference *ref, centre_enum centre):
  d_this(new TranslationalFit_i(ref))
{
  if (centre == fit::cog)
    PositionUtils::translate(&ref->sys(),PositionUtils::cog(ref->sys(),*ref));
  else
    PositionUtils::translate(&ref->sys(),PositionUtils::com(ref->sys(),*ref));
}

TranslationalFit::~TranslationalFit(){
  delete d_this;
}

void TranslationalFit::fit(gcore::System *sys)const{
  PositionUtils::translate(sys,PositionUtils::cog(*sys,*d_this->d_ref));
}

/*
static Vec com(const System &sys, const Reference &ref){
  double totalMass=0;
  Vec cm;

  for(int m=0;m<sys.numMolecules();++m)
    for (int i=0;i < sys.mol(m).numAtoms(); i++) {
      cm += 
	ref.weight(m,i)
	* sys.mol(m).topology().atom(i).mass()
	* sys.mol(m).pos(i) ;

      totalMass += ref.weight(m,i) *
	sys.mol(m).topology().atom(i).mass();

    }

  cm = (1.0/totalMass)*cm;
  return cm;
}


static Vec cog(const System &sys, const Reference &ref){
  Vec cg;

  for(int m=0;m<sys.numMolecules();++m)
    for (int i=0;i < sys.mol(m).numAtoms(); i++) {

     cg += 
	ref.weight(m,i)
	* sys.mol(m).pos(i) ;
    }

  return cg;
}
*/
