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

// pb_FFTBoundaryCondition.cc
#include "FFTBoundaryCondition.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <string>

using pb::FFTBoundaryCondition;

FFTBoundaryCondition::FFTBoundaryCondition(int type, std::string stype,
					   double alpha1, double alpha2, int nalias1, int nalias2, double cutoff, double epsRF, ofstream &os):ppp(os){
  // from public top
  this->tinynum = ppp.getTiny_real();
    
  this->type=type;
  this->stype=stype;
  this->alpha1=alpha1;
  this->alpha2=alpha2;
  this->nalias1=nalias1;
  this->nalias2=nalias2;
  this->cutoff=cutoff;

  this->eps=epsRF;

}



void FFTBoundaryCondition::dumpparameters(ofstream &os) {
  os << "# BOUNDARY PARAMETERS" << endl;
  os << "# -------------------" << endl;
  os << "# TYPE " << type << endl;
  os << "# STYPE " << stype << endl;
  os << "# HAT CHARGE SHAPING FUNCTION" << endl;
  os << "# ALPHA1: " <<  alpha1 << endl;
  os << "# ALPHA2: " << alpha2 << endl;
  os << "# NALIAS1: " << nalias1 << endl;
  os << "# NALIAS2: " << nalias2 << endl;
}
