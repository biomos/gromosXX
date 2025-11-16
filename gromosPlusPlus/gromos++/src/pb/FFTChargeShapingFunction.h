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

// pb_FFTChargeShapingFunction.h

#ifndef INCLUDED_PB_FFTChargeShapingFunction
#define INCLUDED_PB_FFTChargeShapingFunction

#include "PB_Parameters.h"

namespace pb{


class FFTChargeShapingFunction{
public:
  pb::PB_Parameters ppp;	
	int HatType;
	int ParabolaType;
//	int OtherParabolaType;
     
        double tinynum;



public:
    //constructor
      FFTChargeShapingFunction(ofstream &os);
    //deconstructor
     virtual ~FFTChargeShapingFunction(){}
    //methods
	/* FFTChargeShapingFunction getCSF(int type);*/
	virtual double calc(double xx, double yy, double zz,
		            double alpha, double eps0) {return 0.0;}




}; // class
} // namespace


#endif

