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

// pb_FDPoissonIterator.cc

#include "../../config.h"
#ifdef HAVE_LIBFFTW3
#include "FFTPoissonIterator.h"

#include <cmath>
#include <fftw3.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <vector>

#include "FFTGridType.h"
#include "FFTInteractionTypeCodes.h"
#include "FFTBoundaryCondition.h"
#include "FFTDipoleDipole.h"
#include "FFTChargeDipole.h"
#include "FFTDipoleDipole_RF.h"
#include "FFTChargeDipole_RF.h"
#include "FFTDipoleDipole_LS.h"
#include "FFTChargeDipole_LS.h"
#include "PB_Parameters.h"
#include "../fit/PositionUtils.h"
#include "../utils/AtomSpecifier.h"
#include "../gcore/System.h"
#include "../gromos/Exception.h"


using pb::FFTPoissonIterator;
using pb::FFTInteractionTypeCodes;
using pb::FFTChargeDipole;
using pb::FFTDipoleDipole;

FFTPoissonIterator::FFTPoissonIterator(utils::AtomSpecifier atoms, utils::AtomSpecifier atoms_to_charge, int maxsteps, double convergence, double lambda,
				       FFTGridType gt, FFTBoundaryCondition bc, double epssolvent, bool split_potentialbool, ofstream &os):ppp(epssolvent, os), bc(os), gt(os){

  this->atoms = atoms;
  this->atoms_to_charge = atoms_to_charge;
  this->gt = gt;
  this->bc = bc;
  this->epssolvent = epssolvent;
  this->tinynum=ppp.getTiny_real();
  this->convergence=convergence;
  this->maxsteps=maxsteps;
  this->lambda=lambda;
  this->split_potentialbool=split_potentialbool;
                

  this->ls_boundary_eps = ppp.get_ls_boundary_eps();
  this->rf_boundary_eps = ppp.get_rf_boundary_eps();
  this->sc_boundary_eps = ppp.get_sc_boundary_eps();

           
	
	
  FFTInteractionTypeCodes interx_codes;
        
  os << "# FFTPoissonIterator: lambda " << lambda << endl;
  os << "# FFTPoissonIterator: maxsteps " << maxsteps << endl;
  os << "# FFTPoissonIterator: ls_boundary_eps " << ls_boundary_eps << endl;
  os << "# FFTPoissonIterator: rf_boundary_eps " << rf_boundary_eps << endl;
  os << "# FFTPoissonIterator: sc_boundary_eps " << sc_boundary_eps << endl;
  os << "# FFTPoissonIterator: epssolvent " << epssolvent << endl;


  try{
    if (bc.type == interx_codes.lsType) {
      os << "# FFTPoissonIterator: get InterxTensor for LS " << endl;
      ddTensor = new FFTDipoleDipole_LS(ls_boundary_eps,epssolvent, os);
      cdTensor = new FFTChargeDipole_LS(epssolvent, os);
    } else if (bc.type == interx_codes.rfType) {
      os << "# FFTPoissonIterator: get InterxTensor for RF " << endl;
      ddTensor = new FFTDipoleDipole_RF(rf_boundary_eps, bc.cutoff,epssolvent, os);
      cdTensor = new FFTChargeDipole_RF(rf_boundary_eps, bc.cutoff,epssolvent, os);
    } else if (bc.type == interx_codes.scType) {
      os << "# FFTPoissonIterator: get InterxTensor for SC " << endl;
      ddTensor = new FFTDipoleDipole_RF(sc_boundary_eps, bc.cutoff,epssolvent, os);
      cdTensor = new FFTChargeDipole_RF(sc_boundary_eps, bc.cutoff,epssolvent, os);
    }
    else{
      throw gromos::Exception("FFTPoissonIterator","Boundary type not known ...");
    } // end of if
  } // end of try
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }



  os << "# FFTPoissonIterator setup ... make fftplans " << endl;
  setFFTPlans(ppp.get_threadnum(), gt.ngrdx, gt.ngrdy, gt.ngrdz);

}
	
	
	
int FFTPoissonIterator::iterate_poisson(
					std::vector<double> & Vx, std::vector<double> & Vy,  std::vector<double> &Vz,
					std::vector<double> & Ex, std::vector<double> & Ey,  std::vector<double> & Ez,
					std::vector<double> & inside,
					std::vector<double> & RPot, ofstream &os, vector <double> *potentials) {


  int nx=gt.ngrdx; int nx_2=gt.ngrdx/2;
  int ny=gt.ngrdy; int ny_2=gt.ngrdy/2;
  int nz=gt.ngrdz; int nz_2=gt.ngrdz/2;
		
		
  //build k-vectors, store them so we don't need to recalculate during iteration
  std::vector<double> k_vecX;
  std::vector<double> k_vecY;
  std::vector<double> k_vecZ;

  k_vecX.resize(nx, 0.0);
  k_vecY.resize(ny, 0.0);
  k_vecZ.resize(nz, 0.0);


  for (int i=0;i<nx;i++) k_vecX[i] = (((i+nx_2)%nx)-nx_2) * gt.dkx;
  for (int j=0;j<ny;j++) k_vecY[j] = (((j+ny_2)%ny)-ny_2) * gt.dky;
  for (int k=0;k<nz;k++) k_vecZ[k] = (((k+nz_2)%nz)-nz_2) * gt.dkz;
		
  int steps = 1;
  bool converged = false;
  double eFieldInside = 0;
  double prevEFieldInside = 0.0;

  os << "# START ITERATION ... "  <<  endl;
  while (true) { // iteration loop
			
    os << "# ITER "  <<  steps << endl;
    os << "# lambda " << lambda << endl;





    // step 2 in the paper
    fourierTransformedVacuumField(Vx, Vy, Vz, Ex, Ey, Ez);
			
    os << "# done FT vacuum field "  << endl;
			
	
    os << "# FT_VFX " << Vx[0] << " " << Vx[(int)floor(nx_2)] << " " << Vx[nx-1] << endl;
    os << "# FT_VFY " << Vy[0] << " " << Vy[(int)floor(ny_2)] << " " << Vy[ny-1] << endl;
    os << "# FT_VFZ " << Vz[0] << " " << Vz[(int)floor(nz_2)] << " " << Vz[nz-1] << endl;
			
			
			
    // step 3 in the paper
    computeEfieldFromVacuumField(
				 Ex, Ey, Ez, 
				 nx, ny, nz, 
				 k_vecX, k_vecY, k_vecZ, os);

    os << "# done compute Efield from vacuum field "  << endl;

    os << "# FT_EFX " << Ex[0] << " " << Ex[(int)floor(nx_2)] << " " << Ex[nx-1] << endl;
    os << "# FT_EFY " << Ey[0] << " " << Ey[(int)floor(ny_2)] << " " << Ey[ny-1] << endl;
    os << "# FT_EFZ " << Ez[0] << " " << Ez[(int)floor(nz_2)] << " " << Ez[nz-1] << endl;

                  
			
    // step 4 in the paper
    realSpaceElectricField(Ex, Ey, Ez);
			
    os << "# RS_EFX " << Ex[0] << " " << Ex[(int)floor(nx_2)] << " " << Ex[nx-1] << endl;
    os << "# RS_EFY " << Ey[0] << " " << Ey[(int)floor(ny_2)] << " " << Ey[ny-1] << endl;
    os << "# RS_EFZ " << Ez[0] << " " << Ez[(int)floor(nz_2)] << " " << Ez[nz-1] << endl;
			
		
			
    // step 5 in paper
    prevEFieldInside = eFieldInside;
    eFieldInside = computeResidualField(Ex,Ey,Ez,inside, os);
    double deltaE;
    if (1 == steps) {
      deltaE = 0.0;
    } else {
      deltaE = prevEFieldInside - eFieldInside;
    }
			
    if (eFieldInside <= convergence ) {
      converged = true;
      break;
    }
			
    if (steps == maxsteps) break;
			
			
    // this is described in the paper, but we're doing adaptive
    // overrelaxation, so we shouldn't need it.
    //  if ( (eFieldInside-prevEFieldInside > 0.0)&&(steps>1) ) break;
			
    // step 6a in the paper
    updateVacuumField(Vx, Vy, Vz, Ex, Ey, Ez, inside, deltaE, os);
			
    os << "# UP_VFX " << Vx[0] << " " << Vx[(int)floor(nx_2)] << " " << Vx[nx-1] << endl;
    os << "# UP_VFY " << Vy[0] << " " << Vy[(int)floor(ny_2)] << " " << Vy[ny-1] << endl;
    os << "# UP_VFZ " << Vz[0] << " " << Vz[(int)floor(nz_2)] << " " << Vz[nz-1] << endl;
			
    steps++;
  }
  postIteration(Vx, Vy, Vz, Ex, Ey, Ez, RPot, nx, ny, nz, k_vecX, k_vecY, k_vecZ, steps, converged, os, potentials);
		
  return steps;
}
	
	
	
void FFTPoissonIterator::postIteration(
				       std::vector<double> & Vx, std::vector<double> & Vy, std::vector<double> & Vz,
				       std::vector<double> & Ex, std::vector<double> & Ey, std::vector<double> &Ez,
				       std::vector<double> & RPot,
				       int nx, int ny, int nz,
				       std::vector<double> & k_vecX,std::vector<double> & k_vecY,std::vector<double> & k_vecZ,
				       int steps, bool converged, ofstream &os, vector <double> *potentials) {
		
  // step 7 in paper
  // we need the fourier transformed electric field here ...
  // we only do this once, so it's ok, i've timed it. v.
  fft_grid(Ex, Ex);
  fft_grid(Ey, Ey);
  fft_grid(Ez, Ez);

		
  os << "# POST_FT_EFX " << Ex[0] << " " << Ex[(int)floor(nx/2.0)] << " " << Ex[nx-1] << endl;
  os << "# POST_FT_EFY " << Ey[0] << " " << Ey[(int)floor(ny/2.0)] << " " << Ey[ny-1] << endl;
  os << "# POST_FT_EFZ " << Ez[0] << " " << Ez[(int)floor(nz/2.0)] << " " << Ez[nz-1] << endl;
		
		
  reactionFieldHat(
		   Ex, Ey, Ez, 
		   RPot, 
		   nx, ny, nz, 
		   k_vecX, k_vecY, k_vecZ);
  bft_grid(Ex);
  bft_grid(Ey);
  bft_grid(Ez);
  // done step 7, now we need the real-space electric field back.
		
  os << "# POST_RS_EFX " << Ex[0] << " " << Ex[(int)floor(nx/2.0)] << " " << Ex[nx-1] << endl;
  os << "# POST_RS_EFY " << Ey[0] << " " << Ey[(int)floor(ny/2.0)] << " " << Ey[ny-1] << endl;
  os << "# POST_RS_EFZ " << Ez[0] << " " << Ez[(int)floor(nz/2.0)] << " " << Ez[nz-1] << endl;
		
		
		
		
  bft_grid(RPot);
  /* Print average Potental in the inside */
  if (converged) { 
    //os << "# Convergence after " << steps << " steps with dG [kJ/mol]: " << free_energy(RPot) << endl;
    os << "# Convergence after " << steps << " steps with dG [kJ/mol]: " << free_energy_restricted(RPot, os, potentials) << endl;
  }
  else {
    //os << "# Terminated after " << steps << " steps with dG [kJ/mol]: " << free_energy(RPot) << endl;
    os << "# Terminated after " << steps << " steps with dG [kJ/mol]: " << free_energy_restricted(RPot, os, potentials) << endl;
  }					
		
		
		
		
  /* Splitting the reaction potential */
  if (split_potentialbool) {
    split_potential(Vx,Vy,Vz,Ex,Ey,Ez, os);
  }
}
	
	
	
/* step 6a in the paper */
void FFTPoissonIterator::updateVacuumField(
					   std::vector<double> & Vx, std::vector<double> &Vy, std::vector<double> & Vz,
					   std::vector<double> & Ex, std::vector<double> & Ey, std::vector<double> & Ez,
					   std::vector<double> & inside,
					   double deltaSigma, ofstream &os) {
		
		
  int size = gt.ngr3 * 2;
		
  os << "# deltaSigma = " << deltaSigma << endl;


  if (deltaSigma >= 0.0) // we're still converging -- make it run
    lambda *= 1.2;
  else 
    lambda *= 0.5; // we overshot, take it easy now
		
  os << "# Continuing with lambda: " << lambda << endl;
		
  int shift;
  shift = 0;
  for (int index=0;index<size;index+=2) {
			
    if (0.0 != inside[index - shift]) {
				
      double lam_fact = - lambda * inside[index-shift];
      Vx[index] += lam_fact * Ex[index];
      Vy[index] += lam_fact * Ey[index];
      Vz[index] += lam_fact * Ez[index];
      Vx[index+1] += lam_fact * Ex[index+1];
      Vy[index+1] += lam_fact * Ey[index+1];
      Vz[index+1] += lam_fact * Ez[index+1];
    }
    ++shift;
  }
}
	
	
	
	
/* Step 4 in the paper */
	
void FFTPoissonIterator::realSpaceElectricField(
						std::vector<double> & Ex,
						std::vector<double> & Ey,
						std::vector<double> & Ez) {
		
  bft_grid(Ex);
  //enforceZeroImaginaryComponent(Ex);
		
  bft_grid(Ey);
  //enforceZeroImaginaryComponent(Ey);
		
  bft_grid(Ez);
  //enforceZeroImaginaryComponent(Ez);
}
	
void FFTPoissonIterator::enforceZeroImaginaryComponent(std::vector<double> & complexVector) {
		
  for (unsigned int ii = 1; ii < complexVector.size(); ii += 2)
    complexVector[ii] = 0.0;
}
	
	
/*  Step 7 in the paper */
	 
void FFTPoissonIterator::reactionFieldHat(
					  std::vector<double> & Ex,
					  std::vector<double> & Ey,
					  std::vector<double> & Ez,
					  std::vector<double> & RPot,
					  int nx, 
					  int ny, 
					  int nz, 
					  std::vector<double> & k_vecX,
					  std::vector<double> & k_vecY,
					  std::vector<double> & k_vecZ) {
		
  double es_1 = epssolvent- 1;
		
  int shift = 0;
  for (int i=0;i<nx;i++) {
    double kx = k_vecX[i];
    double kx2 = kx*kx;
    for (int j=0;j<ny;j++) {	
      double ky = k_vecY[j];
      double ky2 = ky*ky;
      for (int k=0;k<nz;k++) {
	double kz = k_vecZ[k];
	double kz2 = kz*kz;
	int index = k + nz * ( j + ny * i );
	double k2 = kx2+ky2+kz2;
	
	index+=shift;
	
	if(k2> tinynum) {
	  double pp = cdTensor->polarization(k2) * es_1 / k2;
	  
	  RPot[index] = pp * 
	    (kx*Ex[index+1] + ky*Ey[index+1] + kz*Ez[index+1]);
	  RPot[index+1] = -pp * 
	    (kx*Ex[index] + ky*Ey[index] + kz*Ez[index]);
	}
	else {
	  RPot[index] = 0.0;
	  RPot[index+1] = 0.0;
	}
	++shift;
      }
    }
  }
}
	
	
/* Step 3 in the paper */
	 
void FFTPoissonIterator::computeEfieldFromVacuumField(
						      std::vector<double> & Ex,
						      std::vector<double> & Ey,
						      std::vector<double> & Ez,
						      int nx, 
						      int ny, 
						      int nz, 
						      std::vector<double> & k_vecX,
						      std::vector<double> & k_vecY,
						      std::vector<double> & k_vecZ,
						      ofstream &os) {
		
  double tt[3][3];
		
  int shift = 0;
  for (int i=0;i<nx;i++) {  /* Step 3  */
    double kx = k_vecX[i];
    double kx2 = kx*kx;
    for (int j=0;j<ny;j++) {	
      double ky = k_vecY[j];
      double ky2 = ky*ky;
      for (int k=0;k<nz;k++) {
	double kz = k_vecZ[k];
	double kz2 = kz*kz;
	int index = k + nz * ( j + ny * i );
	double k2 = kx2+ky2+kz2;	
					
	// this may look wasteful, but it actually isn't.
	// trust me. v.
	for (int ii = 0; ii < 3; ii++) tt[0][ii] = kx;
	for (int ii = 0; ii < 3; ii++) tt[1][ii] = ky;
	for (int ii = 0; ii < 3; ii++) tt[2][ii] = kz;
	for (int ii = 0; ii < 3; ii++) tt[ii][0] *= kx;
	for (int ii = 0; ii < 3; ii++) tt[ii][1] *= ky;
	for (int ii = 0; ii < 3; ii++) tt[ii][2] *= kz;


	if (ppp.get_debugvar()==1){
	  os << "# computeEfieldFromVacuumField: tt[0][0] " << tt[0][0] << endl;
	  os << "# computeEfieldFromVacuumField: tt[1][0] " << tt[1][0] << endl;
	  os << "# computeEfieldFromVacuumField: tt[2][0] " << tt[2][0] << endl;
	  os << "# computeEfieldFromVacuumField: tt[0][1] " << tt[0][1] << endl;
	  os << "# computeEfieldFromVacuumField: tt[1][1] " << tt[1][1] << endl;
	  os << "# computeEfieldFromVacuumField: tt[2][1] " << tt[2][1] << endl;
	  os << "# computeEfieldFromVacuumField: tt[0][2] " << tt[0][2] << endl;
	  os << "# computeEfieldFromVacuumField: tt[1][2] " << tt[1][2] << endl;
	  os << "# computeEfieldFromVacuumField: tt[2][2] " << tt[2][2] << endl;
	}
	ddTensor->updateTensor(k2, tt, os);
                                        
	if (ppp.get_debugvar()==1){  
	  os << "# computeEfieldFromVacuumField after upd: tt[0][0] " << tt[0][0] << endl;
	  os << "# computeEfieldFromVacuumField after upd: tt[1][0] " << tt[1][0] << endl;
	  os << "# computeEfieldFromVacuumField after upd: tt[2][0] " << tt[2][0] << endl;
	  os << "# computeEfieldFromVacuumField after upd: tt[0][1] " << tt[0][1] << endl;
	  os << "# computeEfieldFromVacuumField after upd: tt[1][1] " << tt[1][1] << endl;
	  os << "# computeEfieldFromVacuumField after upd: tt[2][1] " << tt[2][1] << endl;
	  os << "# computeEfieldFromVacuumField after upd: tt[0][2] " << tt[0][2] << endl;
	  os << "# computeEfieldFromVacuumField after upd: tt[1][2] " << tt[1][2] << endl;
	  os << "# computeEfieldFromVacuumField after upd: tt[2][2] " << tt[2][2] << endl;
	}



	index+=shift;	
					
                              
                                        
	double er[3];
	er[0]=Ex[index];
	er[1]=Ey[index];
	er[2]=Ez[index];
	double ei[3];
	ei[0]=Ex[index + 1];
	ei[1]=Ey[index + 1];
	ei[2]=Ez[index + 1];
                                        
                                        
	if (ppp.get_debugvar()==1){
	  os << "# index " << index << endl;
	  os << "# er-vec: " << er[0] << " " <<  er[1]<< " " <<  er[2] << endl;
	  os << "# ei-vec: " << ei[0] << " " <<  ei[1]<< " " <<  ei[2] << endl;
	}


	// dot products

	/*	Ex[index]= gmath::Vec.dot(tt[0], er);
		Ey[index]= gmath::Vec.dot(tt[1], er);
		Ez[index]= gmath::Vec.dot(tt[2], er);
		Ex[index+1]=gmath::Vec.dot(tt[0], ei);
		Ey[index+1]=gmath::Vec.dot(tt[1], ei);
		Ez[index+1]=gmath::Vec.dot(tt[2], ei); */


	Ex[index]=0;
	Ey[index]=0;
	Ez[index]=0;
	Ex[index+1]=0;
	Ey[index+1]=0;
	Ez[index+1]=0;

	for (int uu=0; uu<3; uu++){
	  Ex[index]+=tt[0][uu]*er[uu];
	  Ey[index]+=tt[1][uu]*er[uu];
	  Ez[index]+=tt[2][uu]*er[uu];
	  Ex[index+1]+=tt[0][uu]*ei[uu];
	  Ey[index+1]+=tt[1][uu]*ei[uu];
	  Ez[index+1]+=tt[2][uu]*ei[uu];

	}

	if (ppp.get_debugvar()==1){
	  os << "# Ex[index] " << Ex[index] << endl;
	  os << "# Ey[index] " << Ey[index] << endl;
	  os << "# Ez[index] " << Ez[index] << endl;
	  os << "# Ex[index+1] " << Ex[index+1] << endl;
	  os << "# Ey[index+1] " << Ey[index+1] << endl;
	  os << "# Ez[index+1] " << Ez[index+1] << endl;
	}


	++shift;
      }
    }
  }
}
	
	
	
	
/* Step 2 in the paper */
	
void FFTPoissonIterator::fourierTransformedVacuumField(
						       std::vector<double> & Vx,
						       std::vector<double> & Vy,
						       std::vector<double> & Vz,
						       std::vector<double> & Ex,
						       std::vector<double> & Ey,
						       std::vector<double> & Ez) {
  fft_grid(Vx,Ex); 
  fft_grid(Vy,Ey);
  fft_grid(Vz,Ez);
}
	
	
/* Step 5 in the paper */
	
double FFTPoissonIterator::computeResidualField(
						std::vector<double> & Ex, std::vector<double> & Ey, std::vector<double> & Ez,
						std::vector<double> & inside, ofstream &os) {
		
  double E_abs, E_ave = 0.0;
  int count = 0;
  int size = gt.ngr3 * 2;
  int shift = 0;
		
  double imaginaryComponent = 0.0;
  double realComponent = 0.0;
		
  for (int ii=0;ii<size;ii+=2) {
			
    // in the range [0, 1]
    double inSoluteFactor = inside[ii - shift];
			

    //os << "# @*** RESFIELD : ii = " << ii << " inSoluteFactor = " << inSoluteFactor << endl;

    if (inSoluteFactor > tinynum){
      E_abs = 
	Ex[ii]*Ex[ii] +
	Ex[ii+1]*Ex[ii+1] +
	Ey[ii]*Ey[ii] +
	Ey[ii+1]*Ey[ii+1] +
	Ez[ii]*Ez[ii] +
	Ez[ii+1]*Ez[ii+1];
				
      E_ave += E_abs * inSoluteFactor;
      count += inSoluteFactor;			
    }
			
    // FOR DEBUGGING :
    imaginaryComponent += 
      Ex[ii+1]*Ex[ii+1] +
      Ey[ii+1]*Ey[ii+1] +
      Ez[ii+1]*Ez[ii+1];
    realComponent += 
      Ex[ii]*Ex[ii] +
      Ey[ii]*Ey[ii] +
      Ez[ii]*Ez[ii];
					
    // END FOR DEBUGGING
    ++shift;
  } 
		
		
  imaginaryComponent = sqrt(imaginaryComponent);
  realComponent = sqrt(realComponent);

  if (ppp.get_debugvar()==1){
    os << "# imaginary component of the electric field (and average): " <<
      imaginaryComponent << " " << imaginaryComponent/inside.size() << endl;
    os << "# real component of the electric field (and average): " <<
      realComponent << " " << realComponent/inside.size() << endl;
  }
		
  return sqrt(E_ave/count);
}
	
	
	
/********** get the free energy by interpolating *******
 ******* the potential at the location of the charges *******/
double FFTPoissonIterator::free_energy(std::vector<double> & pot) {
		
  double dg;
  int i,j,k;
  unsigned int ion;
  double p;
  double fx,fy,fz;
		
                
		
		
  dg = 0.0; 
  for (ion=0;ion<atoms.size();ion++) {
    // this is just a periodicity-shifting if things are outside the box
    // (should not be necessary in principle if we are using gathering)
    double shiftX = floor( (atoms.pos(ion))[0]/gt.xlen)*gt.xlen;
    // the ion position in terms of cells
    double iPosCX = ( (atoms.pos(ion))[0]-shiftX)/gt.drx;
    // rounded down to the next integer
    i = (int) (iPosCX);
    // offset of the particle with respect to the origin of the cell it's in
    fx = iPosCX - i;
			
    // the same for the y and z coordinates:
    j = (int) ((   (atoms.pos(ion))[1]-floor(  (atoms.pos(ion))[1]/gt.ylen)*gt.ylen)/gt.dry);
    fy =       (   (atoms.pos(ion))[1]-floor(  (atoms.pos(ion))[1]/gt.ylen)*gt.ylen)/gt.dry - j;
    k = (int) ((   (atoms.pos(ion))[2]-floor(  (atoms.pos(ion))[2]/gt.zlen)*gt.zlen)/gt.drz);
    fz =       (   (atoms.pos(ion))[2]-floor(  (atoms.pos(ion))[2]/gt.zlen)*gt.zlen)/gt.drz - k;
			
		
			
			
    // the index of the cell the ion is in (in the pot array)
    int shift0 = k  + gt.ngrdz*(j  +gt.ngrdy*i);
    // the surrounding cells in the positive directions (shifted by 
    // one along one or more of x, y, z).
    int shift1 = k+1+gt.ngrdz*(j  +gt.ngrdy*i    );
    int shift2 = k  +gt.ngrdz*(j+1+gt.ngrdy*i    );
    int shift3 = k+1+gt.ngrdz*(j+1+gt.ngrdy*i    );
    int shift4 = k  +gt.ngrdz*(j  +gt.ngrdy*(i+1));
    int shift5 = k+1+gt.ngrdz*(j  +gt.ngrdy*(i+1)); 
    int shift6 = k  +gt.ngrdz*(j+1+gt.ngrdy*(i+1));
    int shift7 = k+1+gt.ngrdz*(j+1+gt.ngrdy*(i+1));
			
    //System.out.println(i + " " + j + " " + k);
			
    // +     fx *     fy *     fz *pot[k+1+gt.ngrdz*(j+1+gt.ngrdy*(i+1)) + shift7];	
			
    // we multiply the shifts by two, because pot represents a complex vector
    p = 
      (1.0-fx)*(1.0-fy)*(1.0-fz)*pot[shift0 + shift0]
      +(1.0-fx)*(1.0-fy)*     fz *pot[shift1 + shift1]
      +(1.0-fx)*     fy *(1.0-fz)*pot[shift2 + shift2]
      +(1.0-fx)*     fy *     fz *pot[shift3 + shift3]
      +     fx *(1.0-fy)*(1.0-fz)*pot[shift4 + shift4]
      +     fx *(1.0-fy)*     fz *pot[shift5 + shift5]
      +     fx *     fy *(1.0-fz)*pot[shift6 + shift6]
      +     fx *     fy *     fz *pot[shift7 + shift7];
			
			
    dg += atoms.charge(ion) * p;
                        
  }
  return 0.5*dg;
}




/********** get the free energy by interpolating *******
 ******* the potential at the location of the charges *******/
double FFTPoissonIterator::free_energy_restricted(std::vector<double> & pot, ofstream &os, vector<double> *potentials) {


  double dg;
  int i,j,k;
  unsigned int ion;
  double p;
  double fx,fy,fz;




  dg = 0.0;
  for (ion=0;ion<atoms_to_charge.size();ion++) {

    double shiftX = floor( (atoms_to_charge.pos(ion))[0]/gt.xlen)*gt.xlen;
    // the ion position in terms of cells:
    double iPosCX = ( (atoms_to_charge.pos(ion))[0]-shiftX)/gt.drx;
    // rounded down to the next integer
    i = (int) (iPosCX);
    // offset of the particle with respect to the origin of the cell it's in
    fx = iPosCX - i;

    // the same for the y and z coordinates:
    j = (int) ((   (atoms_to_charge.pos(ion))[1]-floor(  (atoms_to_charge.pos(ion))[1]/gt.ylen)*gt.ylen)/gt.dry);
    fy =       (   (atoms_to_charge.pos(ion))[1]-floor(  (atoms_to_charge.pos(ion))[1]/gt.ylen)*gt.ylen)/gt.dry - j;
    k = (int) ((   (atoms_to_charge.pos(ion))[2]-floor(  (atoms_to_charge.pos(ion))[2]/gt.zlen)*gt.zlen)/gt.drz);
    fz =       (   (atoms_to_charge.pos(ion))[2]-floor(  (atoms_to_charge.pos(ion))[2]/gt.zlen)*gt.zlen)/gt.drz - k;




    // the index of the cell the ion is in (in the pot array)
    int shift0 = k  + gt.ngrdz*(j  +gt.ngrdy*i);
    // the surrounding cells in the positive directions (shifted by
    // one along one or more of x, y, z).
    int shift1 = k+1+gt.ngrdz*(j  +gt.ngrdy*i    );
    int shift2 = k  +gt.ngrdz*(j+1+gt.ngrdy*i    );
    int shift3 = k+1+gt.ngrdz*(j+1+gt.ngrdy*i    );
    int shift4 = k  +gt.ngrdz*(j  +gt.ngrdy*(i+1));
    int shift5 = k+1+gt.ngrdz*(j  +gt.ngrdy*(i+1));
    int shift6 = k  +gt.ngrdz*(j+1+gt.ngrdy*(i+1));
    int shift7 = k+1+gt.ngrdz*(j+1+gt.ngrdy*(i+1));

    //System.out.println(i + " " + j + " " + k);

    //+     fx *     fy *     fz *pot[k+1+gt.ngrdz*(j+1+gt.ngrdy*(i+1)) + shift7];

    // we multiply the shifts by two, because pot represents a complex vector
    p =
      (1.0-fx)*(1.0-fy)*(1.0-fz)*pot[shift0 + shift0]
      +(1.0-fx)*(1.0-fy)*     fz *pot[shift1 + shift1]
      +(1.0-fx)*     fy *(1.0-fz)*pot[shift2 + shift2]
      +(1.0-fx)*     fy *     fz *pot[shift3 + shift3]
      +     fx *(1.0-fy)*(1.0-fz)*pot[shift4 + shift4]
      +     fx *(1.0-fy)*     fz *pot[shift5 + shift5]
      +     fx *     fy *(1.0-fz)*pot[shift6 + shift6]
      +     fx *     fy *     fz *pot[shift7 + shift7];



    // print out the electrostatic potential at that atom site
    double pottmp= p;
    if (potentials != NULL) {
      potentials->push_back(pottmp);
    }
    os << "# atom  " << ion << " charge " << atoms_to_charge.charge(ion)  << " ele. pot. " << pottmp << endl;



    dg += atoms_to_charge.charge(ion) * p;

  }
  return 0.5*dg;
}











	
/*
  Splitting the reaction potential into various contributions
  Vx, Vy, Vz are handed over to recompute the last E and 
  as containers for three complex grids for the contributions 
  to the potential 
*/
void FFTPoissonIterator::split_potential(
					 std::vector<double> & Vx, std::vector<double> & Vy, std::vector<double> & Vz,
					 std::vector<double> & Ex, std::vector<double> & Ey, std::vector<double> & Ez, ofstream &os) {
		
  int i,j,k,index;
  int k_vecx,k_vecy,k_vecz;
  double kx2,ky2,kz2,kxy;
  double kx,ky,kz,k2,kr;
  double Tkxx,Tkxy,Tkxz,Tkyy,Tkyz,Tkzz;

  // Maria, here in [0] I store real and in [1] the imag
  double E_tempx[2];
  double E_tempy[2];
  double E_tempz[2];


                
		
  double es_1_es,es_1;
  double c1 = 0.0, c2 = 0.0, Z, cRF = 0.0, cTC = 0.0;  
		
  int nx=gt.ngrdx,nx_2=gt.ngrdx/2;
  int ny=gt.ngrdy,ny_2=gt.ngrdy/2;
  int nz=gt.ngrdz,nz_2=gt.ngrdz/2;
		
		
  if (epssolvent > 0.0) /* "normal" case -> finite epss */
    es_1_es = (epssolvent - 1.0)/epssolvent;
  else es_1_es = 1.0;  /* infinite epss */
  es_1 = epssolvent - 1.0;
		
  fourierTransformedVacuumField(Vx, Vy, Vz, Ex, Ey, Ez);
		
  for (i=0;i<nx;i++) {
    k_vecx = ((i+nx_2)%nx)-nx_2;
    kx = gt.dkx * k_vecx;
    kx2 = kx*kx;
    for (j=0;j<ny;j++) {	
      k_vecy = ((j+ny_2)%ny)-ny_2;
      ky = gt.dky * k_vecy;
      ky2 = ky*ky;
      kxy = kx*ky;
      for (k=0;k<nz;k++) {
	k_vecz = ((k+nz_2)%nz)-nz_2;
	kz = gt.dkz * k_vecz;
	kz2 = kz*kz;
	index = k + nz * ( j + ny * i );
	k2 = kx2+ky2+kz2;
	if (k2<tinynum) {
	  switch (bc.type) {
	  case 2: /* cutoff BC */
	    /* vacuum boundary cond */
	    Tkxx = Tkyy = Tkzz = 1.0 - es_1/(es_1+3.0);
	    Tkxy = Tkyz = Tkxz = 0.0;
	    break;
	  default: /* Barker & Watts RF and EWBC -> both conducting */
	    Tkxx = Tkyy = Tkzz = 1.0; 
	    Tkxy = Tkyz = Tkxz = 0.0;
	    break;
	  }
	} else {
	  /* prefactors c1 and c2 to determine the BC dependent tensor T */
	  switch (bc.type) {
	  case 1: /* EWBC */
	    c1 = 1.0;
	    c2 = - es_1_es;
	    break;
	  case 2: /* cutoff BC */
	    kr = sqrt(k2)*bc.cutoff;
	    Z = 1.0/(kr*kr) * cos(kr) - 1.0/(kr*kr*kr) * sin(kr);	    
	    c1 = 1.0 / (1.0 - es_1 * Z);
	    c2 = - es_1 * (1.0 + 3.0 * Z) / ((1.0 - es_1 * Z) * (1.0 + es_1 * (1.0 + 2.0 * Z)));
	    break;
	  case 3: /* Barker & Watts RF */
	    kr = sqrt(k2)*bc.cutoff;
	    Z = 1 + 3.0/(kr*kr) * cos(kr) - 3.0/(kr*kr*kr) * sin(kr);	
	    c1 = 1.0;
	    c2= - es_1 * Z / (1.0 + es_1 * Z);
	    break;
	  }
	  Tkxx = c1 + c2 * kx2/k2;
	  Tkxy =      c2 * kxy/k2;
	  Tkxz =      c2 * kx*kz/k2;
	  Tkyy = c1 + c2 * ky2/k2;
	  Tkyz =      c2 * ky*kz/k2;
	  Tkzz = c1 + c2 * kz2/k2;
	}
	/*	E_tempx.re[0]=Tkxx*Ex[index]+Tkxy*Ey[index]+Tkxz*Ez[index];
		E_tempx.im[0]=Tkxx*Ex[index+1]+Tkxy*Ey[index+1]+Tkxz*Ez[index+1];
		E_tempy.re[0]=Tkxy*Ex[index]+Tkyy*Ey[index]+Tkyz*Ez[index];
		E_tempy.im[0]=Tkxy*Ex[index+1]+Tkyy*Ey[index+1]+Tkyz*Ez[index+1];
		E_tempz.re[0]=Tkxz*Ex[index]+Tkyz*Ey[index]+Tkzz*Ez[index];
		E_tempz.im[0]=Tkxz*Ex[index+1]+Tkyz*Ey[index+1]+Tkzz*Ez[index+1];
		Ex[index]=E_tempx.re[0];
		Ey[index]=E_tempy.re[0];
		Ez[index]=E_tempz.re[0];
		Ex[index+1]=E_tempx.im[0];
		Ey[index+1]=E_tempy.im[0];
		Ez[index+1]=E_tempz.im[0];*/

	E_tempx[0] =Tkxx*Ex[index]+Tkxy*Ey[index]+Tkxz*Ez[index];
	E_tempx[1]=Tkxx*Ex[index+1]+Tkxy*Ey[index+1]+Tkxz*Ez[index+1];
	E_tempy[0]=Tkxy*Ex[index]+Tkyy*Ey[index]+Tkyz*Ez[index];
	E_tempy[1]=Tkxy*Ex[index+1]+Tkyy*Ey[index+1]+Tkyz*Ez[index+1];
	E_tempz[0]=Tkxz*Ex[index]+Tkyz*Ey[index]+Tkzz*Ez[index];
	E_tempz[1]=Tkxz*Ex[index+1]+Tkyz*Ey[index+1]+Tkzz*Ez[index+1];
	Ex[index]=E_tempx[0];
	Ey[index]=E_tempy[0];
	Ez[index]=E_tempz[0];
	Ex[index+1]=E_tempx[1];
	Ey[index+1]=E_tempy[1];
	Ez[index+1]=E_tempz[1];

	
	/* Reaction potential in k-space via Polarization */	
	switch (bc.type) {
	case 1: /* EWBC */
	  c1 = 1.0;
	  break;
	case 2: /* cutoff BC */
	  kr = sqrt(k2)*bc.cutoff;
	  c1 = 1.0 - sin(kr) / kr;
	  break;  
	case 3: /* BWRF */
	  kr = sqrt(k2)*bc.cutoff;
	  c1 = 1 + 3.0/(kr*kr) * cos(kr) - 3.0/(kr*kr*kr) * sin(kr);
	  cTC = 1.0 - sin(kr) / kr;
	  cRF = c1 - cTC;	
	  break;
	}
	if(k2>tinynum) {
	  Vx[index] = (es_1/k2) * c1 * (kx*Ex[index+1] + ky*Ey[index+1] + kz*Ez[index+1]);
	  Vx[index+1] = (-es_1/k2) * c1 * (kx*Ex[index] + ky*Ey[index] + kz*Ez[index]);
	  Vy[index] = (es_1/k2) * cTC * (kx*Ex[index+1] + ky*Ey[index+1] + kz*Ez[index+1]);
	  Vy[index+1] = (-es_1/k2) * cTC * (kx*Ex[index] + ky*Ey[index] + kz*Ez[index]);
	  Vz[index] = (es_1/k2) * cRF * (kx*Ex[index+1] + ky*Ey[index+1] + kz*Ez[index+1]);
	  Vz[index+1] = (-es_1/k2) * cRF * (kx*Ex[index] + ky*Ey[index] + kz*Ez[index]);
	  //Vincent: looks weird; real component assigned to imaginary and the other way around
	  //shouldn't it be rather:
	  //Vx[index+1] = (es_1/k2) * c1 * (kx*Ex[index+1] + ky*Ey[index+1] + kz*Ez[index+1]);
	  //Vx[index] = (-es_1/k2) * c1 * (kx*Ex[index] + ky*Ey[index] + kz*Ez[index]);
	  //Vy[index+1] = (es_1/k2) * cTC * (kx*Ex[index+1] + ky*Ey[index+1] + kz*Ez[index+1]);
	  //Vy[index] = (-es_1/k2) * cTC * (kx*Ex[index] + ky*Ey[index] + kz*Ez[index]);
	  //Vz[index+1] = (es_1/k2) * cRF * (kx*Ex[index+1] + ky*Ey[index+1] + kz*Ez[index+1]);
	  //Vz[index] = (-es_1/k2) * cRF * (kx*Ex[index] + ky*Ey[index] + kz*Ez[index]);
	}
	else {
	  Vx[index] = 0.0;
	  Vx[index+1] = 0.0;
	  Vy[index] = 0.0;
	  Vy[index+1] = 0.0;
	  Vz[index] = 0.0;
	  Vz[index+1] = 0.0;
	}	  
      }
    }
  }
		
  realSpaceElectricField(Vx, Vy, Vz);
  os << "# contributions to dGsolv (tot, TC, RF): " << free_energy(Vx) << " " << free_energy(Vy) << " " << free_energy(Vz) << endl;
		
  return;
}
	
	
/* Forward direction is needed out of place */
void FFTPoissonIterator::fft_grid(std::vector<double> & r_space, std::vector<double> & k_space) {
  double fact;
  int i;
		
  fact = gt.vol/gt.ngr3;	

  //        os << "# fft_grid: gt.ngr3 " << gt.ngr3 << endl;

  // Maria, copy to a local fftw_complex array

  fftw_complex *in, *out;
  in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*gt.ngr3);
  out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*gt.ngr3);

  //  os << "# fft_grid: assigned in and out fftw vecs " << endl;

  for (int i=0;i<gt.ngr3;i++){
    //      os << "# fft_grid: fill data  -- now i " << i << endl;
    int ii=2*i;
    in[i][0]=r_space[ii];
    in[i][1]=r_space[ii+1];
    out[i][0]=k_space[ii];
    out[i][1]=k_space[ii+1];
  }
  //      os << "# fft_grid: done filling data " << endl;

  //	 my_planV3_f.doFFT(my_planV3_f.plan,r_space,k_space);
  //     fftw_execute(my_planV3_f,r_space,k_space);
  fftw_execute_dft(my_planV3_f,in,out);
  //      os << "# fft_grid: done fftw" << endl;

                
  for (int i=0;i<gt.ngr3;i++){
    //    os << "# fft_grid: fetch data -- now i " << i << endl;
    int ii=2*i;
    r_space[ii]=in[i][0];
    r_space[ii+1]=in[i][1];
    k_space[ii]=out[i][0];
    k_space[ii+1]=out[i][1];
  }
  //     os << "# fft_grid: done fetching data " << endl;

  // free the space
  fftw_free(in);
  fftw_free(out);



  int size = gt.ngr3 * 2;
  for (i=0;i < size;i+=2) {
    k_space[i]   *= fact;
    k_space[i+1] *= fact;
  }		
}
	
/* Backward direction is needed in place */
void FFTPoissonIterator::bft_grid( std::vector<double> & rk_space) {

  
  int i;
		
  double fact = 1.0/gt.vol;

  // Maria, copy to a local fftw_complex array

  fftw_complex *in;
  in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*gt.ngr3);
               


  for (int i=0;i<gt.ngr3;i++){
    int ii=2*i;
    in[i][0]=rk_space[ii];
    in[i][1]=rk_space[ii+1];
       
  }

  // my_planV3_br.doFFT(my_planV3_br.plan,rk_space, rk_space);
  //    fftw_execute(my_planV3_br,rk_space,rk_space);

  fftw_execute_dft(my_planV3_br,in,in);

  for (int i=0;i<gt.ngr3;i++){
    int ii=2*i;
    rk_space[ii]=in[i][0];
    rk_space[ii+1]=in[i][1];
  }


  // free the space
  fftw_free(in);
   


  int size = gt.ngr3 * 2;
  for (i=0;i< size;i+=2) {
    rk_space[i]   *= fact;
    rk_space[i+1] *= fact;
  }
}



void FFTPoissonIterator::setFFTPlans(
				     int numFFTwThreads,
				     int nx, int ny, int nz) {


  // use fftw version 3.x
  //
  // we adopt the array layout as required by FFTW, i.e.
  // element (i) is REAL, wheras
  // element (i+1) is the corresponding COMPLEX to (i).
  // this makes the iteration sometimes a bit cumbersome...

  //std::vector<double>  test1;
  //  std::vector<double>  test2;

  // test1.resize(gt.ngr3 * 2);
  // test2.resize(gt.ngr3 * 2);



  // Maria, use fftw_complex;
  fftw_complex *test1, *test2;
  test1=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*gt.ngr3);
  test2=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*gt.ngr3);
  //    fftw_complex test1[gt.ngr3];
  //   fftw_complex test2[gt.ngr3];



  //	my_planV3_f = new FFTW3DPlanV3(
  //			nx, ny, nz,
  //			test1, test2,
  //			FFTW3DPlanV3.FFTW_FORWARD,
  //			FFTW3DPlanV3.FFTW_MEASURE|FFTW3DPlanV3.FFTW_OUT_OF_PLACE);

  //	my_planV3_br = new FFTW3DPlanV3(
  //			nx, ny, nz,
  //			test1, test2,
  //      		FFTW3DPlanV3.FFTW_BACKWARD,
  //			FFTW3DPlanV3.FFTW_MEASURE|FFTW3DPlanV3.FFTW_OUT_OF_PLACE);

  //  fftw_plan my_planV3_f;
  // fftw_plan my_planV3_br;


  my_planV3_f = fftw_plan_dft_3d(nx,ny,nz,test1,test2,FFTW_FORWARD,FFTW_MEASURE);
  // Maria does it not have to be in place if we reuse the plan??
  // my_planV3_br = fftw_plan_dft_3d(nx,ny,nz,test1,test2,FFT_BACKWARD,FFTW_MEASURE);
  my_planV3_br = fftw_plan_dft_3d(nx,ny,nz,test1,test1,FFTW_BACKWARD,FFTW_MEASURE);

  // free the space
  fftw_free(test1);
  fftw_free(test2);
}

#endif
