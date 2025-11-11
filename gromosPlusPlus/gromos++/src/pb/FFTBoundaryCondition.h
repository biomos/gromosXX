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

// pb_FFTBoundaryCondition.h
#ifndef INCLUDED_PB_FFTBoundaryCondition
#define INCLUDED_PB_FFTBoundaryCondition

#include "PB_Parameters.h"

namespace pb{



class FFTBoundaryCondition{

public:

  pb::PB_Parameters ppp;


        double tinynum;
        
	int type;
	
	std::string stype; 
	/* 
	  Ewald parameter alpha1
	  width of the charge-shaping function
	 */
	double alpha1;
	/*
	  Ewald parameter alpha2
	  width of the charge-shaping function
	 */
	double alpha2;


	/*
	  Ewald parameter charge_shape
	  1 : hat function
          2 : parabola function
          3 : other parabola function
	 
	int charge_shape;
         */

	/*
	  Ewald parameter nalias1
	  number of alias vectors for vacuum potential
	 */

	int nalias1;
	/* 
	  Ewald parameter nalias2
	  number of alias vectors for vacuum potential
	 */
	int nalias2;     
	/*
	  Ewald parameter rd_field
	  controls reading k-space potntial from file
	 
	boolean rd_field = false;*/

	/*
	  cutoff radius
	 */
	double cutoff;          


 /* epsilon (for checking) */
        double eps;

 public:
  // constructor

	FFTBoundaryCondition(int type, std::string stype, double alpha1, double alpha2, int nalias1, int nalias2, double cutoff , double eps, ofstream &os );
 FFTBoundaryCondition(ofstream &os):ppp(os){}

   // deconstructor
    ~FFTBoundaryCondition(){}


  // methods


	void dumpparameters(ofstream &os); 

}; // class
} // namespace


#endif
