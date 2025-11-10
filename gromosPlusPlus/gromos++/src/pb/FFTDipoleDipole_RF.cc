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

// pb_FFTDipoleDipole_RF.cc
#include "FFTDipoleDipole_RF.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>

#include "FFTDipoleDipole_LS.h"
#include "FFTDipoleDipole.h"

using pb::FFTDipoleDipole_RF;
using pb::FFTDipoleDipole;

FFTDipoleDipole_RF::FFTDipoleDipole_RF(
				       double epsilonBoundary,
				       double cutoff, double epsS, ofstream &os): FFTDipoleDipole(epsS,epsilonBoundary, os), ew(epsilonBoundary,epsS, os){
	
  /*    FFTDipoleDipole_RF::FFTDipoleDipole_RF(
	double epsilonBoundary,
	double cutoff){*/

	
  this->epsB=epsilonBoundary;
  this->twoEpsPrimePlusOne = 2 * epsB + 1;
  this->epsSminusOne = epsS - 1;
  this->upper = epsSminusOne * twoEpsPrimePlusOne * twoEpsPrimePlusOne;
  this->twoEs1ePrime = 2 * (epsS - 1) * epsB;
  this->esTwoEprimePlusOne = epsS * twoEpsPrimePlusOne;
  this->cut = cutoff;

  //    if (ppp.get_debugvar()==1){
  os << "# FFTDipoleDipole_RF : epsB " << epsB << endl;
  os << "# FFTDipoleDipole_RF : epsS " << epsS << endl;
  //  }


  // Vincent: we keep one of those around, because
  // their expressions for the k0-vector are the same

  //ew = new pb::FFTDipoleDipole_LS(epsB);
}
	
	
double FFTDipoleDipole_RF::fkr(double kr) {
		
  double kri = 1 / kr;
  double kr2i = kri * kri;
		
  return 3 * kr2i * ( - cos(kr) + sin(kr) * kri);
}
	
double FFTDipoleDipole_RF::computeAFactor(double k2) {
  fckr = fkr(sqrt(k2) * cut);
  if (0.0 == epsB) {
    // conducting boundary
    return 1.0;
  } else {
    return twoEpsPrimePlusOne / (twoEpsPrimePlusOne + epsSminusOne * fckr);
  }
}
	
double FFTDipoleDipole_RF::computeBFactor(double k2) {
  // this assumes computeAFactor(k2) has already been called
  // and fckr is up-to-date. ugly, i agree.
		
  if (0.0 == epsB) {
    // conducting boundary
    return - epsSminusOne * (1 - fckr) / (1 + epsSminusOne * (1 - fckr));
  }
  else {
    return - upper * (fckr - 1) / (
				   (twoEpsPrimePlusOne + epsSminusOne * fckr) *
				   (twoEs1ePrime * fckr - esTwoEprimePlusOne)
				   );
  }
}
	
void FFTDipoleDipole_RF::updateTensorK0(  double  (& tensor)[3][3])  {
  ew.updateTensorK0(tensor);			
}






