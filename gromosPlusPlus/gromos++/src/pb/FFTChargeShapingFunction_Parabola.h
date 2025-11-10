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

// pb_FFTChargeShapingFunction_Parabola.h

#ifndef INCLUDED_PB_FFTChargeShapingFunction_Parabola
#define INCLUDED_PB_FFTChargeShapingFunction_Parabola

#include "FFTChargeShapingFunction.h"

namespace pb{



class FFTChargeShapingFunction_Parabola: virtual public FFTChargeShapingFunction{


public:
        //constructor

	FFTChargeShapingFunction_Parabola(ofstream &os);


        //deconstructor

        ~FFTChargeShapingFunction_Parabola(){};


        // methods

	double calc(double x, double y, double z,
				double alpha, double eps0);



}; // class
} // namespace


#endif
