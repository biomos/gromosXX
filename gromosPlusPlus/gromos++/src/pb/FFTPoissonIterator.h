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

// pb_FFTPoissonIterator.h

#ifndef INCLUDED_PB_FFTPoissonIterator
#define INCLUDED_PB_FFTPoissonIterator
#include <fftw3.h>

#include "FFTGridType.h"
#include "FFTDipoleDipole.h"
#include "FFTChargeDipole.h"
#include "FFTBoundaryCondition.h"
#include "PB_Parameters.h"
#include "../utils/AtomSpecifier.h"

namespace pb{



  class FFTPoissonIterator{
	
	
    utils::AtomSpecifier atoms;
    utils::AtomSpecifier atoms_to_charge;

    pb::PB_Parameters ppp;
    pb::FFTBoundaryCondition bc;
    pb::FFTGridType gt;
    //pb::FFTDipoleDipole *ddTensor;
    //pb::FFTChargeDipole *cdTensor;
        

    double tinynum;
    double convergence;
    int maxsteps;
    double lambda;
    double epssolvent;
    double ls_boundary_eps;
    double rf_boundary_eps;
    double sc_boundary_eps;

    bool split_potentialbool;

    fftw_plan my_planV3_br;
    fftw_plan my_planV3_f;

    FFTDipoleDipole *ddTensor;
    FFTChargeDipole *cdTensor;



  public:
    //constructor
    FFTPoissonIterator(utils::AtomSpecifier atoms, utils::AtomSpecifier atoms_to_charge, int maxsteps, double convergence, double lambda,
		       FFTGridType gt, FFTBoundaryCondition bc, double epssolvent, bool split_potential, ofstream &os);



    // deconstructor
    ~FFTPoissonIterator(){}

    //methods
    void setFFTPlans(
		     int numFFTwThreads,
		     int nx, int ny, int nz);


    int iterate_poisson(
			std::vector<double> & Vx, std::vector<double> & Vy,  std::vector<double> &Vz,
			std::vector<double> & Ex, std::vector<double> & Ey,  std::vector<double> & Ez,
			std::vector<double> & inside,
			std::vector<double> &  RPot,
			ofstream &os, vector <double> *potentials=NULL);


    void postIteration(
		       std::vector<double> & Vx, std::vector<double> & Vy, std::vector<double> & Vz,
		       std::vector<double> & Ex, std::vector<double> & Ey, std::vector<double> &Ez,
		       std::vector<double> & RPot,
		       int nx, int ny, int nz,
		       std::vector<double> & k_vecX,std::vector<double> & k_vecY,std::vector<double> & k_vecZ,
		       int steps, bool converged, ofstream &os, vector <double> *potentials=NULL);


    void updateVacuumField(
			   std::vector<double> & Vx, std::vector<double> &Vy, std::vector<double> & Vz,
			   std::vector<double> & Ex, std::vector<double> & Ey, std::vector<double> & Ez,
			   std::vector<double> & inside,
			   double deltaSigma, ofstream &os);


    void realSpaceElectricField(
				std::vector<double> & Ex,
				std::vector<double> & Ey,
				std::vector<double> & Ez);

    void enforceZeroImaginaryComponent(std::vector<double> & complexVector);

    void reactionFieldHat(
			  std::vector<double> & Ex,
			  std::vector<double> & Ey,
			  std::vector<double> & Ez,
			  std::vector<double> & RPot,
			  int nx,
			  int ny,
			  int nz,
			  std::vector<double> & k_vecX,
			  std::vector<double> & k_vecY,
			  std::vector<double> & k_vecZ);
        
    void computeEfieldFromVacuumField(
				      std::vector<double> & Ex,
				      std::vector<double> & Ey,
				      std::vector<double> & Ez,
				      int nx,
				      int ny,
				      int nz,
				      std::vector<double> & k_vecX,
				      std::vector<double> & k_vecY,
				      std::vector<double> & k_vecZ,
				      ofstream &os);

    void fourierTransformedVacuumField(
				       std::vector<double> & Vx,
				       std::vector<double> & Vy,
				       std::vector<double> & Vz,
				       std::vector<double> & Ex,
				       std::vector<double> & Ey,
				       std::vector<double> & Ez);


    double computeResidualField(
				std::vector<double> & Ex, std::vector<double> & Ey, std::vector<double> & Ez,
				std::vector<double> & inside, ofstream &os);


    double free_energy(std::vector<double> & pot);
    double free_energy_restricted(std::vector<double> & pot, ofstream &os, vector <double> *potentials=NULL);


    void split_potential(
			 std::vector<double> & Vx, std::vector<double> & Vy, std::vector<double> & Vz,
			 std::vector<double> & Ex, std::vector<double> & Ey, std::vector<double> & Ez, ofstream &os);


    void fft_grid(std::vector<double> & r_space, std::vector<double> & k_space);

    void bft_grid(std::vector<double> & rk_space);

  }; // class
} // namespace


#endif
