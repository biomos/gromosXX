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

// pb_FFTChargeShapingFunction.cc
#include "FFTChargeShapingFunction.h"

#include <cstdlib>
#include <cassert>
#include <fstream>

#include "PB_Parameters.h"

using pb::FFTChargeShapingFunction;

FFTChargeShapingFunction::FFTChargeShapingFunction(ofstream &os):ppp(os){

	this->tinynum = ppp.getTiny_real();
      
        this->HatType=0;
        this->ParabolaType=1;
    //    this->OtherParabolaType=2;

 

}

     /*  FFTChargeShapingFunction FFTChargeShapingFunction::getCSF(int type){
       try{
       if (type == HatType){return FFTChargeShapingFunction_Hat();}
       else if (type == ParabolaType){
           return FFTChargeShapingFunction_Parabola();
       }
       else{
         throw gromos::Exception("FFTChargeShapingFunction","Invalid charge shaping function code. Exiting.");
       }
       }// end of try
        catch (const gromos::Exception &e){
                                         cerr << e.what() << endl;
                                         exit(1);
                                }
       }
	*/

	
	
