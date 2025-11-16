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

// pb_FFTGridType.cc
#include "FFTGridType.h"

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cassert>

#include "PB_Parameters.h"


using pb::FFTGridType;

FFTGridType::FFTGridType(int ngrdx, int ngrdy, int ngrdz, double xlen, double ylen, double zlen, int ncubes, ofstream &os):ppp(os) {
 
     // grid point numbers
         this->ngrdx = ngrdx;
	 this->ngrdy = ngrdy;
	 this->ngrdz = ngrdz;
	 this->ngrdxy = ngrdx * ngrdy;
	 this->ngr3 = ngrdx * ngrdy * ngrdz;
      // box edges
	 this->xlen = xlen;
	 this->ylen = ylen;
	 this->zlen = zlen;
      // grid center
         this->centerx=0.5*xlen;
         this->centery=0.5*ylen;
         this->centerz=0.5*zlen;


	 double pi2 = 2 * ppp.getPI();

         // volume of box
	 this->vol=xlen*ylen*zlen;

         /* spacing between r-space grid points */
	 this->drx=xlen/ngrdx;
	 this->dry=ylen/ngrdy;
	 this->drz=zlen/ngrdz;
	 /* spacing between k-space grid points */
	 this->dkx=pi2/xlen;
	 this->dky=pi2/ylen;
	 this->dkz=pi2/zlen;

 	 this->ncubes = ncubes;

	}
void FFTGridType::dumpparameters (ofstream &os) {
	std::cout << "# GRID PARAMETERS" << endl;
	std::cout << "# ---------------" << endl;
        std::cout << "# DIM (GPX GPY GPZ): " << ngrdx << " " << ngrdy << " " << ngrdz << endl;
        std::cout << "# CELL: " << xlen << " " << ylen << " " << zlen << endl;
        std::cout << "# SPACE(R): " << drx << " " << dry << " " << drz << endl;
        std::cout << "# SPACE(K): " << dkx << " " << dky << " " << dkz << endl;
        std::cout << "# CUBES: " << ncubes << endl;
        std::cout << "# VOLUME: " << vol << " nm^3" << endl;

	}
