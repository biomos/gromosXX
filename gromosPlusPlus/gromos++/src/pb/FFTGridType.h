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

// pb_FFTGridType.h

#ifndef INCLUDED_PB_FFTGridType
#define INCLUDED_PB_FFTGridType

#include "PB_Parameters.h"

namespace pb{



class FFTGridType{

public:
  pb::PB_Parameters ppp; // moved back to public
	/*
	  grid dimensions X 
	 */
	int ngrdx;
	/* 
	  grid dimensions Y 
	 */
	int ngrdy;
	/* 
	  grid dimensions Z 
	 */
	int ngrdz;
	/* 
	  grid dimensions X * grid dimensions Y 
	 */
	int ngrdxy;	
	/*
	  grid dimensions X * grid dimensions Y * grid dimensions Z
	 */
	int ngr3;
	/*
	  X edge of the periodic unit cell
	 */
	double xlen;
	/* 
	  Y edge of the periodic unit cell
	 */
	 double ylen;
	/*
	  Z edge of the periodic unit cell
	 */
	double zlen;
	/*
	  unit cell volume
	 */
	double vol;


        /* the grid center: the middle of the grid*/
        double centerx;
        double centery;
        double centerz;



	/* 
	  X spacing between r-space grid points
	 */               
	 double drx;
	/*
	  Y spacing between r-space grid points
	 */
	double dry;
	/*
	  Z spacing between r-space grid points
	 */
	double drz;
	/*
	  X spacing between k-space grid points
	 */
	double dkx;
	/*
	  Y spacing between k-space grid points
	 */
	double dky;
	/* 
	  Z spacing between k-space grid points
	 */
	double dkz;
	/*
	  number of cubes
	 */
	int ncubes;
	
	
	 //recycles Ewald vacuum field
	
	//boll recycleVACfield = false;
	
	  //writes Ewald vacuum fields to file
	
	//bool writeVACfield = false;
	
	 // reads Ewald vacuum fields to file
	 
	// boolreadVACfield = false;
	// string vacfield1 = "";
	// string vacfield2 = "";
	
	//bool havevacfield = false;



	



    //constructor

	
	FFTGridType(int ngrdx, int ngrdy, int ngrdz, double xlen, double ylen, double zlen, int ncubes, ofstream &os);
 FFTGridType(ofstream &os):ppp(os){}
    // deconstructor
       ~FFTGridType(){}

        //methods

	void dumpparameters(ofstream &os);

        
}; // class
} // namespace


#endif
