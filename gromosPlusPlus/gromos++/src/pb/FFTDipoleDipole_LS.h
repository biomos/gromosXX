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

// pb_FFTDipoleDipole_LS.h

#ifndef INCLUDED_PB_FFTDipoleDipole_LS
#define INCLUDED_PB_FFTDipoleDipole_LS

#include "FFTDipoleDipole.h"

namespace pb{


class FFTDipoleDipole_LS: public FFTDipoleDipole{

	
	double k0EpsilonFactor;
	double es_1_es;
        double epsB;
        
public:
        //constructor

	FFTDipoleDipole_LS(double epsB,double epsS, ofstream &os);

        //deconstructor

        ~FFTDipoleDipole_LS(){};


        // methods
        
	void updateTensorK0(double (& tensor)[3][3]);
	
	double computeAFactor(double k2);
	
	double computeBFactor(double k2);


}; // class
} // namespace


#endif
