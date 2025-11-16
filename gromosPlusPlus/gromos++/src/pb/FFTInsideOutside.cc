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

// pb_FFTInsideOutside.cc
#include "FFTInsideOutside.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <vector>

#include "PB_Parameters.h"
#include "../fit/PositionUtils.h"
#include "../utils/AtomSpecifier.h"
#include "../gcore/System.h"

using pb::FFTInsideOutside;

FFTInsideOutside::FFTInsideOutside(int ngridpoints, ofstream &os){
  pb::PB_Parameters ppp(os);
  this->tinynum=ppp.getTiny_real();
  this->ngridpoints=ngridpoints;
}



void FFTInsideOutside::inside_sol(
				  int  gridN[3],  double  gridD[3],
				  int nc,
				  utils::AtomSpecifier & atoms,
				  std::vector<double> & in) {
		
		
  // to keep track of grid points on the solute surface
  // cell diameter / 2
  double dr = 0.5 * sqrt(gridD[0]*gridD[0]+gridD[1]*gridD[1]+gridD[2]*gridD[2]);

  // stuff relating to cubelet dimensions
  int ce = 2 * nc + 1;
  double cei = 1.0 / ce;
   
  double ddx = gridD[0] * cei;
  double ddy = gridD[1] * cei;
  double ddz = gridD[2] * cei;
  double ctoi = cei * cei * cei;
		
  // atomic data
  int numAtoms = atoms.size();
	
		
  // loop over grid cells
  for (int ii = 0; ii < gridN[0]; ii++) {

    for (int jj = 0; jj < gridN[1]; jj++) {

      for (int kk = 0; kk < gridN[2]; kk++) {
	int index = kk + gridN[2] * (jj + gridN[1] * ii);

	// for this gridpoint, make an array[(ce)^3] which saves (bool true or false)
	// whether the cubelet is inside an atom
	int cecece = ce*ce*ce;
	bool cear[cecece];
	// initialise it to false
	for (int uu=0;uu<cecece; uu++ ){
	  cear[uu]= false;
	}


	// for this gridpoint, initialise the cin-counter (counts how many cublets are inside an atom) to 0
	int cin = 0;


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




	// loop over all atoms
	for (int aa = 0; aa < numAtoms; aa++) {

                  

	  // skip atom with zero radius
	  double rr = atoms.radius(aa);

	  // go to next atom if no radius
	  if (fabs(rr)<tinynum){
	    continue;
	  }
                        

	  double rplus = rr + dr;
	  double rminus;
	  if (rr > dr) {
	    rminus = rr - dr;
	  }
	  else {
	    rminus = 0.0;
	  }
	  double r2plus = rplus*rplus;
	  double r2minus = rminus*rminus;
	  double r2 = rr*rr;
			
	  double thispos[3];
	  // int thisgridpos[3];
	  thispos[0]=(atoms.pos(aa))[0];
	  thispos[1]=(atoms.pos(aa))[1];
	  thispos[2]=(atoms.pos(aa))[2];

	  // determine the nearest image distance between the atom thispos[] and the grid cell origin
	  // the coordinates of the atom are thispos[]
	  double nimvec[3];
	  nimvec[0]=cellori[0]-thispos[0];
	  nimvec[1]=cellori[1]-thispos[1];
	  nimvec[2]=cellori[2]-thispos[2];

	  for (int u=0; u<3;u++){
	    if (  nimvec[u] > boxh[u]    ) {nimvec[u] = nimvec[u] - box[u];}
	    else if ( nimvec[u]< (-1.0*boxh[u]) ) {nimvec[u] = nimvec[u] + box[u];}
	    else {;}
	  }

	  // compute the distance between the atom and the center
	  // of the other cell
	  double distance2 = nimvec[0]*nimvec[0]+nimvec[1]*nimvec[1]+nimvec[2]*nimvec[2];



								
						
	  if (distance2 >= (r2plus+tinynum)) {
	    // no overlap, go to next atom
	    continue;
	  }


	  // alright, there's going to be _some_ interaction.
	  if ((distance2 < r2minus+ tinynum  )&&(distance2 > tinynum)) {
	    // inside cutoff, count completely
	    cin = cecece;
	    for (int uu=0;uu<cecece; uu++ ){
	      cear[uu]= true;
	    }



	    //	fullyInside++;
	    // print the coordinates of the atom and the grid point
	  }
	  else if (distance2 <= tinynum)  {
	    // atom sits on grid point
	    cin = cecece;
	    for (int uu=0;uu<cecece; uu++ ){
	      cear[uu]= true;
	    }
                                                  

	  }
	  else {
	    int cubeletcounter = 0;

	    for (int xx=-nc;xx<=nc;xx++) {
	      double dx = nimvec[0] + xx * ddx;
	      double dx2 = dx * dx;
								
	      for (int yy=-nc;yy<=nc;yy++) {
		double dy = nimvec[1] + yy * ddy;
		double dy2 = dy * dy;
									
		for (int zz=-nc;zz<=nc;zz++) {
		  double dz = nimvec[2] + zz * ddz;
		  double dz2 = dz * dz;
										
		  double dist2c = dx2 + dy2 + dz2;


                                                                               

		  if (dist2c <  r2 + tinynum) {
											
		    // is this cubelet already taken by another atom?

		    if (cear[cubeletcounter] == false) {
		      // not yet occupied
		      cin++;
		      cear[cubeletcounter] = true;
		    } // end if

		  } // end if
                                                                                
                                                                                
		  cubeletcounter++ ; // end of this cubelet
                                                                                
		}
	      }
	    } // end loop over cubelets
							
                                                        
                                                        
	  }//} // end of else  (cube at cutoff)
	
	} // end of atom-loop

                                   


	double overlap;
	overlap = cin * ctoi; // it is >=1 if all cubelets are occupied
                                  
	in[index] = overlap;
	if ( in[index] > 1.0 ) {
	  in[index] = 1.0;
	}
	// this is necessary because of numerical issues otherwise
	// (for the heaviside function later, when we subtract 1, it should be exactly 1 now so that it gets exactly zero then)
	else if (fabs(in[index] - 1.0) < tinynum ) {
	  in[index] = 1.0;
	}
	else{
	  ;
	}




	// finished with this grid point

      }				
    }
  } // end of grid point loops
				
}
	
	

/* Integrates the inside grid and sets values to 1 or 0,
   computes integral of inside and boundary and returns nr of points inside */
void FFTInsideOutside::integr_inside(std::vector <double>  &  inside, ofstream &os) {



  unsigned int index;
  double integral=0.0;
  int inPoints_noB=inside.size(); /* without boundary */

  /* Eliminate everything from the inside that is not 1.0
     and compute integrals*/
  for (index=0;index< inside.size() ;index++) {
    integral += inside[index];
    if (inside[index] < 1.0-tinynum) {

      // if (  inside[index] > tinynum ){

      inside[index] = 0.0;
      inPoints_noB--;
    }
  }
  os << "# FFTInsideOutside::integr_inside: gridpoints: " <<  inside.size() << endl;
  os << "# FFTInsideOutside::integr_inside: inPoints_noB: " <<  inPoints_noB << endl;
  os << "# FFTInsideOutside::integr_inside: integral " <<     integral << endl;
  os << "# FFTInsideOutside::integr_inside: normalized integral " << integral/(1.0* inside.size()) << endl;

}
