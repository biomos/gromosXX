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

// pb_FFTChargeDipole_LS.h

#ifndef INCLUDED_PB_FFTChargeDipole_LS
#define INCLUDED_PB_FFTChargeDipole_LS

#include "FFTChargeDipole.h"

namespace pb{



class FFTChargeDipole_LS: virtual public FFTChargeDipole{


public:
        //constructor

  FFTChargeDipole_LS(double epssolvent, ofstream &os);


        //deconstructor

        ~FFTChargeDipole_LS(){}


        // methods

	double polarization(double k2);



}; // class
} // namespace


#endif


