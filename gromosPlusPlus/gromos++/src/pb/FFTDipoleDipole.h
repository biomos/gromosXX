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

// pb_FFTDipoleDipole.h

#ifndef INCLUDED_PB_FFTDipoleDipole
#define INCLUDED_PB_FFTDipoleDipole

#include "PB_Parameters.h"


namespace pb{



  class FFTDipoleDipole{

  protected:
    
    double tinynum;
    double epsS;
    double epsB;
    pb::PB_Parameters ppp;

  public:

    //constructor
    FFTDipoleDipole(double epsilonSolvent, double epsilonB, ofstream &os);

    // deconstructor
    virtual ~FFTDipoleDipole(){}

    //methods
    
    
    /* Determine the inverse k-space dipole-dipole interaction
       tensor for the k=0-vector. you shouldn't need this directly*/

    virtual void updateTensorK0(double (& tensor)[3][3]){}


    /* determine the scalar multiplier A of the identity matrix
       component of the interaction tensor. you shouldn't need
       this directly*/

    virtual double computeAFactor(double k2){return 0.0;}

    /* determine the scalar multiplier B of the outer product
       component of the interaction tensor. you shouldn't need
       this directly*/

    virtual double computeBFactor(double k2){return 0.0;}


    /* Determine the inverse k-space dipole-dipole interaction
       tensor*/

    void updateTensor(		double k2,     double  (& tensor)[3][3], ofstream &os
				) ;
    
    
    
    
  }; // class
} // namespace


#endif
