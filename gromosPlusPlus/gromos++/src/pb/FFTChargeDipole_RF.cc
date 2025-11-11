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

// pb_FFTChargeDipole_RF.cc
#include "FFTDipoleDipole_RF.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>

#include "FFTChargeDipole_RF.h"
#include "FFTChargeDipole.h"
#include "../gromos/Exception.h"

using pb::FFTChargeDipole_RF;
using pb::FFTChargeDipole;
using pb::FFTDipoleDipole_RF;

FFTChargeDipole_RF::FFTChargeDipole_RF(
				       double epsilonBoundary,
				       double cutoff,double epssolvent, ofstream &os) : FFTChargeDipole(epssolvent, os) {


  this->epsB=epsilonBoundary;
  this->cut = cutoff;


  os << "# FFTChargeDipole_RF: epssolvent " << epssolvent << endl;
  os << "# FFTChargeDipole_RF: epsB " << epsB << endl;;


  try{
    if (epsB < 1.0 &&   (fabs(epsB) >  tinynum) ) {
      throw gromos::Exception("FFTChargeDipole_RF","Invalid permittivity for boundary (should be 0, 1 or >1) ...");

    }// endif
  }// end of try

                         
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }       
                
                
  if (fabs(epsB)<tinynum) {
    this->first = 1.0;
    this->second = 0.0;
  }
  else {
    double twoEpsPrimePlusOne = 2 * epsB + 1;
    double twoEpsPrimeMinusOne = 2 * (epsB - 1);
			
    this->first = twoEpsPrimeMinusOne / twoEpsPrimePlusOne;
    this->second = 3 / twoEpsPrimePlusOne;
  }			

}

/* kr the product of the norm of the k-vector and the cutoff distance */

double FFTChargeDipole_RF::fkr(double kr) {
  //double res;
  //return FFTDipoleDipole_RF.fkr(kr);


  double kri = 1 / kr;
  double kr2i = kri * kri;

  return 3 * kr2i * ( - cos(kr) + sin(kr) * kri);



}
	
/*  we assume that the A-factor (and therefore fckr) is up-to-date */
	 
double FFTChargeDipole_RF::polarization(double k2) {
		
  double kr = sqrt(k2) * cut;		
		
  return 1.0 - first * fkr(kr) - second * sin(kr) / kr;
}



