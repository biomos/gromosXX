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

// pb_FFTInsideOutside.h

#ifndef INCLUDED_PB_FFTInsideOutside

#include <vector>

#include "PB_Parameters.h"
#include "../utils/AtomSpecifier.h"



namespace pb{



  class FFTInsideOutside{


    double tinynum;
    int ngridpoints;
    //std::vector<double>  in;



  public:
    //constructor
    FFTInsideOutside(int ngridpoints, ofstream &os);



    // deconstructor
    ~FFTInsideOutside(){}

    //methods

    /* determine the grid points in the solute,
       This is identical to inside_sol_orig, but
       the outer loop is the loop over the solute atoms, this way, only grid points
       which may actually be inside the solute need to be investigated, greatly speeding
       up the algorithm */


    void inside_sol(  int gridN[3],  double  gridD[3],
		      int nc,
		      utils::AtomSpecifier & atoms, std::vector<double> & in);


    /* Integrates the inside grid and sets values to 1 or 0,
       computes integral of inside and boundary and returns nr of points inside */
    void integr_inside(std::vector <double>  &  inside, ofstream &os);


  }; // class
} // namespace


#endif
