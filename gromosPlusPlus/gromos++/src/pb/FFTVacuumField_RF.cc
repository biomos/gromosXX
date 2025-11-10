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

// pb_FFTVacuumField_RF.cc
#include "FFTVacuumField_RF.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <vector>

#include "FFTBoundaryCondition.h"
#include "FFTGridType.h"
#include "FFTVacuumField.h"
#include "../gromos/Exception.h"
#include "../fit/PositionUtils.h"
#include "../utils/AtomSpecifier.h"
#include "../gcore/System.h"
//#include "FFTGridUtils.h"


using pb::FFTVacuumField_RF;
using pb::FFTVacuumField;
//using pb::FFTGridUtils;

FFTVacuumField_RF::FFTVacuumField_RF(utils::AtomSpecifier atoms, FFTGridType gt, FFTBoundaryCondition bc, ofstream &os):FFTVacuumField(atoms, gt, bc, os){
  //FFTVacuumField_RF::FFTVacuumField_RF():FFTVacuumField(atoms, gt, bc){

}
	
	
void FFTVacuumField_RF::calcVacField(
				     std::vector<double>& fldx,
				     std::vector<double>& fldy,
				     std::vector<double>& fldz,
				     ofstream &os) {
		
  calcVacFieldVincent(fldx, fldy, fldz, bc.eps, os);
}	
	

void FFTVacuumField_RF::calcVacFieldVincent(
					    std::vector<double> & fldx,
					    std::vector<double> & fldy,
					    std::vector<double> & fldz,
					    double epsRf, ofstream &os) {
		
  os << "# computing vacuum field with epsilon " << epsRf << endl;
  os << "# Warning: The vacuum field on a grid point hosting a particle will be set to zero ..." << endl;
		
  // we need this for the reaction field
  double rcut_3 = 1 / (bc.cutoff*bc.cutoff*bc.cutoff);
  double srf;
  if (fabs(epsRf)<tinynum) // 0: conducting reaction field
    srf = 1;
  else 
    srf = 2 * (epsRf - 1) / (2 * epsRf + 1);
  double FRfs = srf * rcut_3;
		
  // cell diameter / 2
  double dr = 0.5 * sqrt(gt.drx*gt.drx+gt.dry*gt.dry+gt.drz*gt.drz);
		
  double gridLengths[3];
  gridLengths[0]=gt.xlen;
  gridLengths[1]=gt.ylen;
  gridLengths[2]=gt.zlen;

  // cutoff checks
  try{
    if (bc.cutoff < dr){
      throw gromos::Exception("FFTVacuumField_RF","Cutoff too small for this grid (cutoff<dr). Exiting.");
    } 


    if (  (bc.cutoff > gridLengths[0]) || (bc.cutoff > gridLengths[1]) || (bc.cutoff > gridLengths[2])   ) {
      throw gromos::Exception("FFTVacuumField_RF","Cutoff too large for this grid (cutoff>gridlength). Exiting.");
    }
  }// end of try
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
             
		
		
  // the number of grid cells per dimension
  int nc = gt.ncubes;
		
  // grid data
  double gridD[3];
  int gridN[3];

  gridD[0]=gt.drx;
  gridD[1]=gt.dry;
  gridD[2]=gt.drz;
  gridN[0]=gt.ngrdx;
  gridN[1]=gt.ngrdy;
  gridN[2]=gt.ngrdz;
                      
                
  // stuff relating to cubelet dimensions
  double ddx = gt.drx/(2.0*nc+1.0);
  double ddy = gt.dry/(2.0*nc+1.0);
  double ddz = gt.drz/(2.0*nc+1.0);
  int ce = 2 * nc + 1;
  int cto = ce * ce * ce;
  double inverseCubeletsPerCell = 1.0 / cto;

  // cutoff information
  double rcut = bc.cutoff;
  double rplus = rcut + dr;
  double rminus = rcut - dr;
  double r2plus = rplus*rplus;
  double r2minus = rminus*rminus;
  double r2 = rcut*rcut;

  // the cutoff in terms of grid cells
  std::vector <double> scaledCharges;
  scaledCharges.resize(atoms.size());
	
		
  for (unsigned int ii = 0; ii < atoms.size(); ii++) {
    scaledCharges[ii] = fpepsi * atoms.charge(ii);
  }

		
  // set the field to zero everywhere
		
  // NOTE:
  // In the following, the loops over grid points are implemented as loops from 0 to (exclusive) gridN[i],
  // i.e. the loops cover gridN points in each dimension. Note that the length of a grid cell is
  // box_length/gridN, which divides the box edge in gridN grid edges. Counting the grid points,
  // one will see that there are actually gridN+1 points. But: This ist just fine like that. It 
  // has to be like that, because the gridN+first point is already the picture of the zeroeth
  // grid point.
  // The FFT routines work with a periodicity of box_length. In the view of this,
  // the periodicity convention explained above is perfectly fine, too.
		
  for (int ii = 0; ii < gridN[0]; ii++) {
    for (int jj = 0; jj < gridN[1]; jj++) {
      for (int kk = 0; kk < gridN[2]; kk++) {
	//index multiplied by 2; see below for explanation
	int index = 2 * (kk + gridN[2] * (jj + gridN[1] * ii));
	//this calculation only involves REAL numbers...
	//as we store both RE and IM numbers in 1 array
	//(separated by 1 integer, i.e. 
	// num[0].re = num[0]
	// num[0].im = num[1]
	// num[1].re = num[2]
	// num[1].im = num[3]...)
	//this layout is used by FFTW
	//we therefore need to shift up...
	fldx[index]=0.0;fldx[index+1]=0.0;
	fldy[index]=0.0;fldy[index+1]=0.0;
	fldz[index]=0.0;fldz[index+1]=0.0;
      }
    }
  }
		
  // atomic data
  int numAtoms = atoms.size();

  //      FFTGridUtils gridut;

  // loop over all atoms
  for (int aa = 0; aa < numAtoms; aa++) {
			
			
    // skip uncharged atom
    if ( fabs(scaledCharges[aa]) < tinynum) continue;
			
    double pos[3];
    //int gPos[3];
                        
    pos[0]=(atoms.pos(aa))[0];
    pos[1]=(atoms.pos(aa))[1];
    pos[2]=(atoms.pos(aa))[2];

                       
    // loop over grid cells
    for (int ii = 0; ii < gridN[0]; ii++) {

      for (int jj = 0; jj < gridN[1]; jj++) {

	for (int kk = 0; kk < gridN[2]; kk++) {


	  int    index = 2 * (kk + gt.ngrdz * ( jj + gt.ngrdy * ii ));


	  // determine the nearest image distance between the atom pos[] and the grid cell origin



	  // the coordinates of the atom are pos[]
	  // the coordinates of the grid cell origin are gridD * ii, gridD * jj, gridD * kk



	  double cellori[3];
	  cellori[0]=ii*gridD[0];
	  cellori[1]=jj*gridD[1];
	  cellori[2]=kk*gridD[2];

	  double boxh[3];
	  double box[3];
	  boxh[0]=gridN[0]*gridD[0]/2.0;
	  boxh[1]=gridN[1]*gridD[1]/2.0;
	  boxh[2]=gridN[2]*gridD[2]/2.0;
	  box[0]=gridN[0]*gridD[0];
	  box[1]=gridN[1]*gridD[1];
	  box[2]=gridN[2]*gridD[2];

	  double nimvec[3];
	  nimvec[0]=cellori[0]-pos[0];
	  nimvec[1]=cellori[1]-pos[1];
	  nimvec[2]=cellori[2]-pos[2];

	  for (int u=0; u<3;u++){
	    if (  nimvec[u] > boxh[u]    ) {nimvec[u] = nimvec[u] - box[u];}
	    else if ( nimvec[u]< (-1.0*boxh[u]) ) {nimvec[u] = nimvec[u] + box[u];}
	    else {;}
	  }

	  double distance2 = nimvec[0]*nimvec[0]+nimvec[1]*nimvec[1]+nimvec[2]*nimvec[2];
	  double overlap;

                                   

						
	  if (distance2 >= (r2plus+tinynum)  ) {
	    // no overlap, go to next cell
	    continue;
	  }

                    


	  if ( distance2 <= tinynum ){
	    // particle sits on grid point .... set vacuum field to zero
	    fldx[index]=0;
	    fldy[index]=0;
	    fldz[index]=0;
	    continue;
	  }

	  // compute distances 
	  double dist3i = 1.0 / (distance2*sqrt(distance2));
	  double dfact = scaledCharges[aa]*(dist3i - FRfs);
	  double dx,dy,dz;
	  if ((distance2 < r2minus+ tinynum  )&&(distance2 > tinynum)) {
	    /* this cube is fully inside the cutoff */
	    overlap = 1.0;
	  }
	  else { /* this cube is at the cutoff */
	    //   if ((distance2 < r2plus+tinynum)&&(distance2 > tinynum)) {
	    /* run over the little cubes */
	    int    cin = 0;
	    for (int iii=-nc;iii<=nc;iii++) {
	      for (int jjj=-nc;jjj<=nc;jjj++) {
		for (int kkk=-nc;kkk<=nc;kkk++) {
		  dx = nimvec[0] + iii * ddx;
		  dy = nimvec[1]+ jjj * ddy;
		  dz = nimvec[2]+ kkk * ddz;
		  distance2 = dx*dx+dy*dy+dz*dz;
		  if (distance2 < r2+tinynum){
		    cin++;
		  }
		}
	      }
	    }
	    overlap =  cin * inverseCubeletsPerCell;
	  }//}
	  // done with the else if (fully in cube or partly in cube)
                                                                
	  double efact = dfact * overlap;
	  fldx[index] += efact * nimvec[0];
	  fldy[index] += efact * nimvec[1];
	  fldz[index] += efact * nimvec[2];
	
	}				
      }
    } // done, loop over grid
		
  } // done, loop over atoms
       

}
