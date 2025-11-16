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

// pb_FDPoisson.cc

#include "../../config.h"
#ifdef HAVE_LIBFFTW3
#include <fftw3.h>

#include "FFTPoisson.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <iomanip>

#include "../fit/PositionUtils.h"
#include "../utils/AtomSpecifier.h"
#include "../gcore/System.h"

#include "FFTGridType.h"
#include "FFTBoundaryCondition.h"
#include "FFTDipoleDipole.h"
#include "PB_Parameters.h"
#include "FFTInteractionTypeCodes.h"
#include "FFTVacuumField.h"
#include "FFTVacuumField_LS.h"
#include "FFTVacuumField_RF.h"
#include "FFTInsideOutside.h"
#include "../gromos/Exception.h"


using pb::FFTPoisson;
using pb::FFTVacuumField;
using pb::FFTInteractionTypeCodes;
using pb::FFTBoundaryCondition;
using pb::FFTInsideOutside;
using pb::FFTDipoleDipole;

FFTPoisson::FFTPoisson(utils::AtomSpecifier atoms,utils::AtomSpecifier atoms_to_charge, FFTGridType gt, FFTBoundaryCondition bc, double gridspacing, int maxsteps, double convergence, double lambda,
		       double epssolvent, bool split_potentialbool, bool shift_atoms, ofstream &os):
  ppp(epssolvent, os), bc(os),
  pbiterator(atoms, atoms_to_charge,maxsteps, convergence, lambda,  gt, bc, epssolvent,split_potentialbool, os), gt(os)
{

  this->atoms=atoms;
  this->atoms_to_charge=atoms_to_charge;
  this->tinynum=ppp.getTiny_real();
  this->convergence=convergence;
  this->maxsteps=maxsteps;
  this->lambda=lambda;

  this->gt=gt;
  this->bc=bc;

  this->epssolvent=epssolvent;
  this->gridspacing=gridspacing;
  this->split_potentialbool=split_potentialbool;

  // read the radii
  for (unsigned int i=0; i<atoms.size(); i++){
    this->radii.push_back(atoms.radius(i));
  }
  

  os << "# epssolvent " << epssolvent << endl;
  for (unsigned int i=0; i<atoms_to_charge.size();i++){
    os << "# atom " << i << " charge " << atoms_to_charge.charge(i) << endl;
    os << "# atom " << i << " radius " << atoms_to_charge.radius(i) << endl;
  }

		
  //pbiterator = new pb::FFTPoissonIterator(atoms, maxsteps, convergence, lambda,
  //gt, bc, epssolvent, split_potential);
		
  // os << "# FFTPoisson setup ... make fftplans " << endl;
  //setFFTPlans(1, gt.ngrdx, gt.ngrdy, gt.ngrdz);

             
  if (shift_atoms == 1) {
    os << "# FFTPoisson setup ... center atoms " << endl;
    center_atoms_on_grid(atoms, gt.centerx,gt.centery,gt.centerz, os);
  }

}
	

/* void FFTPoisson::setFFTPlans(
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
                                
   fftw_plan my_planV3_f;
   fftw_plan my_planV3_br;


   my_planV3_f = fftw_plan_dft_3d(nx,ny,nz,test1,test2,FFTW_FORWARD,FFTW_MEASURE);
   // Maria does it not have to be in place if we reuse the plan??
   // my_planV3_br = fftw_plan_dft_3d(nx,ny,nz,test1,test2,FFT_BACKWARD,FFTW_MEASURE);
   my_planV3_br = fftw_plan_dft_3d(nx,ny,nz,test1,test1,FFTW_BACKWARD,FFTW_MEASURE);

   // free the space
   fftw_free(test1);
   fftw_free(test2);


   } */
	
void FFTPoisson::solve_poisson(ofstream &os, vector <double> *potentials) {

  //possible scaling of radii of atoms that are close to the border of the box...
  atomshift(os);
  
  //check grid dimension....
  gridcheck(os);
	
		
	     
  int gridN[3];
  double gridD[3];
  gridN[0]=gt.ngrdx;
  gridN[1]=gt.ngrdy;
  gridN[2]=gt.ngrdz;
  gridD[0]=gt.drx;
  gridD[1]=gt.dry;
  gridD[2]=gt.drz;

                
  FFTInsideOutside inandout(gt.ngr3, os);
           
  std::vector<double>  inside;
  inside.resize(gt.ngr3);


  //     for (int iii =0; iii < inside.size() ;iii++) {
  //          inside[iii] = 0.0 ;
  //      }

  os << "# FFTPoisson setup ...  insideoutside" << endl;
  inandout.inside_sol(gridN, gridD, gt.ncubes, atoms, inside);
  // force garbage collection to free surface points hash
  //System.gc();
		
  os << "# FFTPoisson setup ... done insideoutside" << endl;

  //we adopt the array layout as required by FFTW, i.e.
  //element (i) is REAL, wheras
  //element (i+1) is the corresponding COMPLEX to (i).
  //this makes the iteration sometimes a bit cumbersome...
  os << "# FFTPoisson setup ... grid" << endl;
  std::vector<double> Vx;
  std::vector<double> Vy;
  std::vector<double> Vz;

  Vx.resize(gt.ngr3 * 2);
  Vy.resize(gt.ngr3 * 2);
  Vz.resize(gt.ngr3 * 2);

  os << "# FFTPoisson setup ... done gridresizing" << endl;
	     
  os << "# FFTPoisson setup ... setup vacfield " << endl;
  setupVacuumField(inside, Vx, Vy, Vz, os);
  // force garbage collection to free vacuum field recycling stuff if necessary
  //System.gc();
		
		
		
  /* Now eliminate everything from the inside that is not 1.0 */
  inandout.integr_inside(inside, os);
		
  //init grids...
  //we adopt the array layout as required by FFTW, i.e.
  //element (i) is REAL, wheres 
  //element (i+1) is the corresponding COMPLEX to (i).
  //this makes the iteration sometimes a bit cumbersome...


  std::vector<double> Pot;
  std::vector<double> fldx;
  std::vector<double> fldy;
  std::vector<double> fldz;


  Pot.resize(gt.ngr3 * 2);
  fldx.resize(gt.ngr3 * 2);
  fldy.resize(gt.ngr3 * 2);
  fldz.resize(gt.ngr3 * 2);


		
  pbiterator.iterate_poisson(
			     Vx, Vy, Vz,
			     fldx, fldy, fldz, inside, 
			     Pot, os, potentials);
		
}

/* set up initial (unmodified) vacuum field */
	  
/* inside non-zero if the corresponding grid cell is inside the solute
   vx vacuum field - x component
   vy vacuum field - y component
   vz vacuum field - z component
*/


         
void FFTPoisson::setupVacuumField(
				  std::vector<double> & inside,
				  std::vector<double> & vx, std::vector<double> & vy, std::vector<double> & vz, ofstream &os) {

  FFTVacuumField *vacfield;//(atoms, gt, bc);
  FFTInteractionTypeCodes interx_codes;

  os << "# setting up vacuum field ... " << endl;
  os << "# bc.type is: " << bc.type << endl;

  try{
    if (bc.type == interx_codes.lsType) {
      os << "# need LS field " <<  endl;
      vacfield = new FFTVacuumField_LS(atoms,gt,bc, os);
    } else if (bc.type == interx_codes.rfType) {
      os << "# need RF field " <<  endl;
      vacfield = new FFTVacuumField_RF(atoms,gt,bc, os);
    }else if (bc.type == interx_codes.scType) {
      os << "# need SC field" <<  endl;
      vacfield = new FFTVacuumField_RF(atoms,gt,bc, os);
    } else{
      throw gromos::Exception("FFTPoisson","Boundary type not known ...");
    } //end of if
  } // end of try
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  if (   bc.type !=  interx_codes.lsType )
    { 
      vacfield->calcVacField(vx, vy, vz, os);
    }

  else {

    // ewald-type vacuum field
    vacfield->calcVacField(vx, vy, vz, os);
    // we obtain the fourier transformed
    // vacuum field in the ewald case, but we want
    // the r-space vacuum field, so transform:
    pbiterator.bft_grid(vx);
    pbiterator.bft_grid(vy);
    pbiterator.bft_grid(vz);

  }
			






  /* initial guess for the modified vacuum field -
   * set V zero inside of solute */
  int size = gt.ngr3 * 2;
  int shift = 0;
  for (int index=0;index < size;index+=2) {
    double heaviside = (1-inside[index - shift]);



    vx[index] *= heaviside;
    vy[index] *= heaviside;
    vz[index] *= heaviside;
    vx[index+1] *= heaviside;
    vy[index+1] *= heaviside;
    vz[index+1] *= heaviside;


    //  os << "# @@@ vvv " <<  vx[index] << " " <<  vy[index] << " "  << vz[index] << endl;
    //                      os << "# @@@ vvv " <<  vx[index+1] << " " <<  vy[index+1] << " "  << vz[index+1] << endl;


    ++shift;
  }

  os<< "# after setupfield: vx[0] " << vx[0] << endl;
  os<< "# after setupfield: vy[0] " << vy[0] << endl;
  os<< "# after setupfield: vz[0] " << vz[0] << endl;
  os<< "# after setupfield: vx[100] " << vx[100] << endl;
  os<< "# after setupfield: vy[100] " << vy[100] << endl;
  os<< "# after setupfield: vz[100] " << vz[100] << endl;



               
}

/*
void FDPoissonBoltzmann::increasebox(ofstream &os){
  int size = atoms.size();
  double gridposition_X_lower = 0.0;
  double gridposition_X_upper = gt.xlen;
  double gridposition_Y_lower = 0.0;
  double gridposition_Y_upper = gt.ylen;
  double gridposition_Z_lower = 0.0;
  double gridposition_Z_upper = gt.zlen;
  double protruding_X = 0;
  double protruding_Y = 0;
  double protruding_Z = 0;
  double temp_lower = 0;
  double temp_upper = 0;

  for (int i=0; i < size; ++i) {
    // X AXIS
    temp_lower = gridposition_X_lower - atoms.pos(i)[0];
    temp_upper = atoms.pos(i)[0] - gridposition_X_upper;
    if ( temp_lower > protruding_X ) {
      protruding_X = temp_lower;
    }
    if ( temp_upper > protruding_X ) {
      protruding_X = temp_upper;
    }

    // Y AXIS
    temp_lower = gridposition_Y_lower - atoms.pos(i)[1];
    temp_upper = atoms.pos(i)[1] - gridposition_Y_upper;
    if ( temp_lower > protruding_Y ) {
      protruding_Y = temp_lower;
    }
    if ( temp_upper > protruding_Y ) {
      protruding_Y = temp_upper;
    }

        // Z AXIS
    temp_lower = gridposition_Z_lower - atoms.pos(i)[2];
    temp_upper = atoms.pos(i)[2] - gridposition_Z_upper;
    if ( temp_lower > protruding_Z ) {
      protruding_Z = temp_lower;
    }
    if ( temp_upper > protruding_Z ) {
      protruding_Z = temp_upper;
    }
  }
  if (protruding_X < gridspacing) {
    GPX += 1;
    os << "Increased the number of gridpoints in the X direction by 1 because one or more atoms were too close to the border" << endl;
  } else {
    throw gromos::Exception("FDPoissonBolzmann","Atom extending grid. Exiting ...");
  }
  if (protruding_Y < gridspacing) {
    GPY += 1;
    os << "Increased the number of gridpoints in the Y direction by 1 because one or more atoms were too close to the border" << endl;
  } else {
    throw gromos::Exception("FDPoissonBolzmann","Atom extending grid. Exiting ...");
  }
    if (protruding_Z < gridspacing) {
    GPX += 1;
    os << "Increased the number of gridpoints in the Z direction by 1 because one or more atoms were too close to the border" << endl;
  } else {
    throw gromos::Exception("FDPoissonBolzmann","Atom extending grid. Exiting ...");
  }
}
*/

void FFTPoisson::atomshift(ofstream &os) {
  int size = atoms.size();

  for (int i=0; i < size; ++i) {
    os << "# GRIDCHECK for atom " << i+1 << " radius " << fixed << std::setprecision(7) << radii[i]
       << " charge " << fixed << std::setprecision(3) << setw(6) << atoms.charge(i)
       << " pos " << atoms.pos(i)[0] << " " << atoms.pos(i)[1] << " " << atoms.pos(i)[2] << endl;
    
      bool gridfail = 0;
      char coordfail = '0';
      double atomposition = 0;
      double gridposition = 0;
      double protruding = 0;
      
      //X-axis upper end atom exceeding
      gridposition = gt.xlen;
      if ( atoms.pos(i)[0] > gridposition ) {
	coordfail='X';
	protruding = atoms.pos(i)[0] - gridposition;
	atoms.pos(i)[0] = 0.0+protruding;
	os << "# Needed to shift atom " << i+1 << " on the " << coordfail << " axis to the other side of the box " << endl;
	os << "# New coordinates atom " << i+1  << " " << atoms.pos(i)[0] << " " << atoms.pos(i)[1] << " " << atoms.pos(i)[2] << endl;
      };

      //X-axis upper end vdw sphere exceeding
      gridposition = gt.xlen;
      if ( atoms.pos(i)[0] + radii[i] > gridposition) {
	coordfail='X';
	protruding = atoms.pos(i)[0] + radii[i] - gridposition + gridspacing/100;
	radii[i] = radii[i] - protruding;
	os << "# Needed to shrink the radius of atom " << i+1 << " on the " << coordfail << " axis since it exceeded the grid " << endl;
	os << "# New radius atom " << i+1  << " " << radii[i] << endl;
      };

      //X-axis lower end atom exceeding
      gridposition = 0.0;
      if ( atoms.pos(i)[0] < gridposition ) {
	coordfail='X';
	protruding = atoms.pos(i)[0] - gridposition;
	atoms.pos(i)[0] = gt.xlen + protruding;
	os << "# Needed to shift atom " << i+1 << " on the " << coordfail << " axis to the other side of the box " << endl;
	os << "# New coordinates atom " << i+1  << " " << atoms.pos(i)[0] << " " << atoms.pos(i)[1] << " " << atoms.pos(i)[2] << endl;
      };

      //X-axis lower end vdw sphere exceeding
      gridposition = 0.0;
      if ( atoms.pos(i)[0] - radii[i] < gridposition ) {
	coordfail='X';
	protruding = atoms.pos(i)[0] - radii[i] - gridposition - gridspacing/100;
	radii[i] = radii[i] + protruding;
	os << "# Needed to shrink the radius of atom " << i+1 << " on the " << coordfail << " axis since it exceeded the grid " << endl;
	os << "# New radius atom " << i+1  << " " << radii[i] << endl;
      };

      //X-axis upper end vdw sphere exceeding (another one, for the case it was shifted from the lower to the upper end, but radius still exceeds)
      gridposition = gt.xlen;
      if ( atoms.pos(i)[0] + radii[i] > gridposition) {
	coordfail='X';
	protruding = atoms.pos(i)[0] + radii[i] - gridposition + gridspacing/100;
	radii[i] = radii[i] - protruding;
	os << "# Needed to shrink the radius of atom " << i+1 << " on the " << coordfail << " axis since it exceeded the grid " << endl;
	os << "# New radius atom " << i+1  << " " << radii[i] << endl;
      };




      //Y-axis upper end atom exceeding
      gridposition = gt.ylen;
      if ( atoms.pos(i)[1] > gridposition) {
	coordfail='Y';
	protruding = atoms.pos(i)[1] - gridposition;
	atoms.pos(i)[1] = 0.0+protruding;
	os << "# Needed to shift atom " << i+1 << " on the " << coordfail << " axis to the other side of the box " << endl;
	os << "# New coordinates atom " << i+1  << " " << atoms.pos(i)[0] << " " << atoms.pos(i)[1] << " " << atoms.pos(i)[2] << endl;
      };

      //Y-axis upper end vdw sphere exceeding
      gridposition = gt.ylen;
      if ( atoms.pos(i)[1] + radii[i] > gridposition) {
	coordfail='Y';
	protruding = atoms.pos(i)[1] + radii[i] - gridposition + gridspacing/100;
	radii[i] = radii[i] - protruding;
	os << "# Needed to shrink the radius of atom " << i+1 << " on the " << coordfail << " axis since it exceeded the grid " << endl;
	os << "# New radius atom " << i+1  << " " << radii[i] << endl;
      };

      //Y-axis lower end atom exceeding
      gridposition = 0.0;
      if ( atoms.pos(i)[1] < gridposition ) {
	coordfail='Y';
	protruding = atoms.pos(i)[1] - gridposition;
	atoms.pos(i)[1] = gt.ylen + protruding;
	os << "# Needed to shift atom " << i+1 << " on the " << coordfail << " axis to the other side of the box " << endl;
	os << "# New coordinates atom " << i+1  << " " << atoms.pos(i)[0] << " " << atoms.pos(i)[1] << " " << atoms.pos(i)[2] << endl;
      };

      //Y-axis lower end vdw sphere exceeding
      gridposition = 0.0;
      if ( atoms.pos(i)[1] - radii[i] < gridposition ) {
	coordfail='Y';
	protruding = atoms.pos(i)[1] - radii[i] - gridposition - gridspacing/100;
	radii[i] = radii[i] + protruding;
	os << "# Needed to shrink the radius of atom " << i+1 << " on the " << coordfail << " axis since it exceeded the grid " << endl;
	os << "# New radius atom " << i+1  << " " << radii[i] << endl;
      };

      //Y-axis upper end vdw sphere exceeding (another one, for the case it was shifted from the lower to the upper end, but radius still exceeds)
      gridposition = gt.ylen;
      if ( atoms.pos(i)[1] + radii[i] > gridposition ) {
	coordfail='Y';
	protruding = atoms.pos(i)[1] + radii[i] - gridposition + gridspacing/100;
	radii[i] = radii[i] - protruding;
	os << "# Needed to shrink the radius of atom " << i+1 << " on the " << coordfail << " axis since it exceeded the grid " << endl;
	os << "# New radius atom " << i+1  << " " << radii[i] << endl;
      };
      

      //Z-axis upper end atom exceeding
      gridposition = gt.zlen;
      if ( atoms.pos(i)[2] > gridposition ) {
	coordfail='Z';
	protruding = atoms.pos(i)[2] - gridposition;
	atoms.pos(i)[2] = 0.0+protruding;
	os << "# Needed to shift atom " << i+1 << " on the " << coordfail << " axis to the other side of the box " << endl;
	os << "# New coordinates atom " << i+1  << " " << atoms.pos(i)[0] << " " << atoms.pos(i)[1] << " " << atoms.pos(i)[2] << endl;
      };

      //Z-axis upper end vdw sphere exceeding
      gridposition = gt.zlen;
      if ( atoms.pos(i)[2] + radii[i] > gridposition) {
	coordfail='Z';
	protruding = atoms.pos(i)[2] + radii[i] - gridposition + gridspacing/100;
	radii[i] = radii[i] - protruding;
	os << "# Needed to shrink the radius of atom " << i+1 << " on the " << coordfail << " axis since it exceeded the grid " << endl;
	os << "# New radius atom " << i+1  << " " << radii[i] << endl;
      };

      //Z-axis lower end atom exceeding
      gridposition = 0.0;
      if ( atoms.pos(i)[2] < gridposition) {
	coordfail='Z';
	protruding = atoms.pos(i)[2] - gridposition;
	atoms.pos(i)[2] = gt.zlen + protruding;
	os << "# Needed to shift atom " << i+1 << " on the " << coordfail << " axis to the other side of the box " << endl;
	os << "# New coordinates atom " << i+1  << " " << atoms.pos(i)[0] << " " << atoms.pos(i)[1] << " " << atoms.pos(i)[2] << endl;
      };

      //Z-axis lower end vdw sphere exceeding
      gridposition = 0.0;
      if ( atoms.pos(i)[2] - radii[i] < gridposition ) {
	coordfail='Z';
	protruding = atoms.pos(i)[2] - radii[i] - gridposition - gridspacing/100;
	radii[i] = radii[i] + protruding;
	os << "# Needed to shrink the radius of atom " << i+1 << " on the " << coordfail << " axis since it exceeded the grid " << endl;
	os << "# New radius atom " << i+1  << " " << radii[i] << endl;
      };

      //Z-axis upper end vdw sphere exceeding (another one, for the case it was shifted from the lower to the upper end, but radius still exceeds)
      gridposition = gt.zlen;
      if ( atoms.pos(i)[2] + radii[i] > gridposition) {
	coordfail='Z';
	protruding = atoms.pos(i)[2] + radii[i] - gridposition + gridspacing/100;
	radii[i] = radii[i] - protruding;
	os << "# Needed to shrink the radius of atom " << i+1 << " on the " << coordfail << " axis since it exceeded the grid " << endl;
	os << "# New radius atom " << i+1  << " " << radii[i] << endl;
      };
      

      
      
  }// end of for loop
}// end of atomshift

void FFTPoisson::gridcheck(ofstream &os) {
		
  int size = atoms.size();
  double rad;
		
  for (int i=0; i < size; ++i) {
    rad = radii[i];
    try{
      if (   (  (atoms.pos(i))[0] + rad  ) > (gt.xlen)
	     || ((atoms.pos(i))[0] - rad) < (0.0)
	     || ((atoms.pos(i))[1] + rad) > (gt.ylen)
	     || ((atoms.pos(i))[1] - rad) < (0.0)
	     || ((atoms.pos(i))[2]+ rad) > (gt.zlen)
	     || ((atoms.pos(i))[2] - rad) < (0.0) ) {




	os << "# ATOM EXTENDING GRID!" << endl;
	os << "# ATOM ID:     " << (i+1)  << endl;
	os << "# ATOM POS:    " << (atoms.pos(i))[0] << " " << (atoms.pos(i))[1] << " " << (atoms.pos(i))[2]<< endl;
	os << "# ATOM RADIUS: "  <<  rad << endl;
	os << "# GRID START   "  <<  0.0 << " " << 0.0 << " " <<  0.0 << endl;
	os << "# GRID BOX     "  <<  (gt.xlen) << " " << (gt.ylen) << " " << (gt.zlen) << endl;


	throw gromos::Exception("FDPoissonBolzmann","Atom extending grid. Exiting ...");

      } // if endif
    } // end of try

    catch (const gromos::Exception &e){
      cerr << e.what() << endl;
      exit(1);
    }

  } // loop over atoms

}



// shift atoms so that their cog coincides with the grid center
void FFTPoisson::center_atoms_on_grid(utils::AtomSpecifier  & atoms,
				      double gridcenterx, double gridcentery, double gridcenterz, ofstream &os){
//  double cogAtoms[3];
  // get cog
//  cogAtoms[0]=fit::PositionUtils::cog(*atoms.sys(), atoms)[0];
//  cogAtoms[1]=fit::PositionUtils::cog(*atoms.sys(), atoms)[1];
//  cogAtoms[2]=fit::PositionUtils::cog(*atoms.sys(), atoms)[2];

      double mmcenterAtoms[3];
      gmath::Vec coormin = fit::PositionUtils::getmincoordinates(atoms.sys(), false);
      gmath::Vec coormax = fit::PositionUtils::getmaxcoordinates(atoms.sys(), false);
      for (int i=0;  i<3; i++) 
	mmcenterAtoms[i] = (coormin[i] + coormax[i]) * 0.5; 
      

  //os << "# initial cog x,y,z = " <<  cogAtoms[0] << " "<< cogAtoms[1] << " " << cogAtoms[2] << endl;
  os << "# initial min max center x,y,z = " <<  mmcenterAtoms[0] << " "<< mmcenterAtoms[1] << " " << mmcenterAtoms[2] << endl; 
  // get the translation vector
  double tvec[3];
  tvec[0]=gridcenterx-mmcenterAtoms[0];
  tvec[1]=gridcentery-mmcenterAtoms[1];
  tvec[2]=gridcenterz-mmcenterAtoms[2];
  
  os << "translation vector x,y,z = " <<tvec[0] << " " << tvec[1] << " " << tvec[2] << endl;

  // now shift


  for (unsigned int i=0;i<atoms.size(); i++){
    atoms.pos(i)[0] += tvec[0];
    atoms.pos(i)[1] += tvec[1];
    atoms.pos(i)[2] += tvec[2];

  }

  os << "# gridcenter x,y,z = " << gridcenterx << " " << gridcentery << " " << gridcenterz << endl;
  // check the new cog
//  cogAtoms[0]=fit::PositionUtils::cog(*atoms.sys(), atoms)[0];
//  cogAtoms[1]=fit::PositionUtils::cog(*atoms.sys(), atoms)[1];
//  cogAtoms[2]=fit::PositionUtils::cog(*atoms.sys(), atoms)[2];
//  os << "# final cog (after centering on grid) x,y,z = " <<  cogAtoms[0] << " "<< cogAtoms[1] << " " << cogAtoms[2] << endl;
// check the new mmcenter
      coormin = fit::PositionUtils::getmincoordinates(atoms.sys(), false);
      coormax = fit::PositionUtils::getmaxcoordinates(atoms.sys(), false);
      for (int i=0;  i<3; i++)
          mmcenterAtoms[i] = (coormin[i] + coormax[i]) * 0.5;
      os << "# final mmcenterAtoms (after centering on grid) x,y,z = " <<  mmcenterAtoms[0] << " "<< mmcenterAtoms[1] << " " << mmcenterAtoms[2] << endl;

}
#endif
