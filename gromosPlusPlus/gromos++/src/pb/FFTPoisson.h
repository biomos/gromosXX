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

// pb_FFTPoisson.h

#ifndef INCLUDED_PB_FFTPoisson
#define INCLUDED_PB_FFTPoisson

#include "PB_Parameters.h"
#include "FFTPoissonIterator.h"
#include "FFTBoundaryCondition.h"
#include "FFTGridType.h"
#include "FFTInteractionTypeCodes.h"
#include "FFTVacuumField.h"
#include "../utils/AtomSpecifier.h"

namespace pb{



  class FFTPoisson{
	
    utils::AtomSpecifier atoms;
    utils::AtomSpecifier atoms_to_charge;
        
    pb::PB_Parameters ppp;
    pb::FFTPoissonIterator pbiterator;
    pb::FFTBoundaryCondition bc;
    pb::FFTGridType gt;
	
    double tinynum;
    double convergence;
    int maxsteps;
    double lambda;
    double epssolvent;
    double gridspacing;
    bool split_potentialbool;
    vector<double> radii;
    //static j3DFFT j3DFFT;
	
	
	
    //plans for FFTW(V3)
    //fftw_plan my_planV3_f; //forward
    //fftw_plan my_planV3_br; //backward
  public:
    // constructor
    FFTPoisson(utils::AtomSpecifier atoms,utils::AtomSpecifier atoms_to_charge, FFTGridType gt, FFTBoundaryCondition bc, double gridspacing, int maxsteps, double convergence, double lambda,
	       double epssolvent, bool split_potentialbool, bool shift_atoms, ofstream &os);



    // deconstructor
    ~FFTPoisson(){}

    // methods
    // void setFFTPlans(
    //		int numFFTwThreads,
    //			int nx, int ny, int nz);

	
    void solve_poisson(ofstream &os, vector <double> *potentials=NULL);
		
    void setupVacuumField(
			  std::vector<double> &  inside,
			  std::vector<double> & vx,std::vector<double> &  vy, std::vector<double> &  vz, ofstream &os);

    void center_atoms_on_grid(utils::AtomSpecifier  & atoms, double gridcenterx, double gridcentery, double gridcenterz, ofstream &os);

    void atomshift(ofstream &os);
    void gridcheck(ofstream &os);

  }; // class
} // namespace


#endif
