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

// pb_FFTChargeShaping_Function_Hat.cc
#include "FFTChargeShapingFunction_Hat.h"

#include <cmath>
#include <cstdlib>
#include <cassert>
#include <fstream>

#include "FFTChargeShapingFunction.h"

using pb::FFTChargeShapingFunction_Hat;
using pb::FFTChargeShapingFunction;

FFTChargeShapingFunction_Hat::FFTChargeShapingFunction_Hat(ofstream &os):FFTChargeShapingFunction(os){
// FFTChargeShapingFunction_Hat::FFTChargeShapingFunction_Hat(){
}

	

		
		 /* k-space version of the potential due to
		 a spherical hat function at the origin plus homogeneous background charge 
		 Phi(k) = 1/(eps0*k^2) * 12(2 - 2 Math.cos ak - ak Math.sin ak)/(ak)^4 */
		 
	double FFTChargeShapingFunction_Hat::calc(
				double x, double y, double z, 
				double alpha, double eps0) {
			
			
		double k2 = x*x + y*y + z*z;
			if (k2 >= tinynum) {
			double alphak = alpha*sqrt(k2);
				return 1.0/(eps0*k2)*12.0*(2.0-2.0*cos(alphak)-alphak*sin(alphak))/(alphak*alphak*alphak*alphak);
			} else {
				return 0.0;
			}
		}
	


