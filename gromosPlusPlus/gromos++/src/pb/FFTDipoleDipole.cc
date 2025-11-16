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

// pb_FFTDipoleDipole.cc
#include "FFTDipoleDipole.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>

#include "../gromos/Exception.h"
#include "PB_Parameters.h"




using pb::FFTDipoleDipole;

FFTDipoleDipole::FFTDipoleDipole(double epsilonSolvent, double epsilonB, ofstream &os):ppp(epsilonSolvent, os){
	
  this->tinynum = ppp.getTiny_real();
  this->epsS = epsilonSolvent;
  this->epsB=epsilonB;
                

  os << "# FFTDipoleDipole: epsS " << epsS << endl;
  os << "# FFTDipoleDipole: epsB " << epsB << endl;;

                 
  try{
    if (  (epsB < 1.0 && fabs(epsB) >  tinynum)   )  {
      throw gromos::Exception("FFTDipoleDipole","Invalid permittivity for boundary (should be 0, 1 or >1) ...");
    }// endif
  }// end of try

  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }

}
	

	
/* Determine the inverse k-space dipole-dipole interaction
   tensor*/
	
void FFTDipoleDipole::updateTensor(
				   double k2,
				   double (& tensor)[3][3],
				   ofstream &os
				   ) {


  // os << "# FFTDipoleDipole::updateTensor ... " << endl;
		
  if (k2 < tinynum) {
    // this is the zero-vector
    updateTensorK0(tensor);
    //    os << "# was zero vec ..." << endl;
    return;
  }
		
  double A = computeAFactor(k2);
		
  double B = computeBFactor(k2) / k2;
  if (ppp.get_debugvar()==1){
    os << "# A and B " << A << " " << B << endl;
  }

  for (int ii = 0; ii < 3; ii++){
    for (int jj = 0; jj < 3; jj++){
      tensor[ii][jj] *= B;
      if (ppp.get_debugvar()==1){
	os << "# Bfac ... tensor [ " << ii << " ][ " << jj << "] = " << tensor[ii][jj] << endl;
      }
    }}
		
  for (int ii = 0; ii < 3; ii++){
    tensor[ii][ii] += A;
    if (ppp.get_debugvar()==1){
      os << "# Afac ... tensor [ " << ii << " ][ " << ii << "] = " << tensor[ii][ii] << endl;
    }
  }
                
                
}

