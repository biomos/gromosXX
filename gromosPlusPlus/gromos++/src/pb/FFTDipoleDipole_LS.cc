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

// pb_FFTDipoleDipole_LS.cc
#include "FFTDipoleDipole_LS.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cassert>

#include "../gromos/Exception.h"
#include "FFTDipoleDipole.h"



using pb::FFTDipoleDipole_LS;
using pb::FFTDipoleDipole;

FFTDipoleDipole_LS::FFTDipoleDipole_LS(double epsBoundary,double epsS, ofstream &os): FFTDipoleDipole(epsS,epsBoundary, os){
//FFTDipoleDipole_LS::FFTDipoleDipole_LS(double epsBoundary){

		this->es_1_es = (epsS - 1) / epsS;

                this->epsB=epsBoundary;

  try{
		if (epsB < 1.0 &&  (fabs(epsB) >  tinynum)  ) {
			  throw gromos::Exception("FFTDipoleDipole",
                                  "Invalid permittivity for boundary (should be 0, 1 or <1) ...");

		                }// endif
}// end of try
 catch (const gromos::Exception &e){
                                         cerr << e.what() << endl;
                                         exit(1);
                                }
  

		if (fabs(epsB) < tinynum) {
			// conducting boundary
			k0EpsilonFactor = 1.0;
		} else {
			k0EpsilonFactor = (2 * epsB + 1) / (2 * epsB + epsS);
		}




}



	void FFTDipoleDipole_LS::updateTensorK0(double ( & tensor) [3][3]) {
		
		for (int ii = 0; ii < 3; ii++) 
			for (int jj = 0; jj < 3; jj++)
				tensor[ii][jj] = 0.0;
		
		for (int ii = 0; ii < 3; ii++){
			tensor[ii][ii] = k0EpsilonFactor;
                
                }
		
	}
	
	double FFTDipoleDipole_LS::computeAFactor(double k2) {
		return 1.0;
	}
	
	double FFTDipoleDipole_LS::computeBFactor(double k2) {
		return - es_1_es;
	}
