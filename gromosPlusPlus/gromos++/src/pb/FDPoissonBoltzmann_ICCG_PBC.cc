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

// pb_FDPoissonBoltzmann_ICCG_PBC.cc
#include "FDPoissonBoltzmann_ICCG_PBC.h"

#include <cmath>
#include <cstdlib>
#include <cassert>
#include <vector>

using pb::FDPoissonBoltzmann_ICCG_PBC;

FDPoissonBoltzmann_ICCG_PBC::FDPoissonBoltzmann_ICCG_PBC(int GPX, int GPY, int GPZ){



		this->GPX = GPX;
		this->GPY = GPY;
		this->GPZ = GPZ;
		this->GPXGPY = GPX*GPY;
		this->GPXGPYGPZ = GPX*GPY*GPZ;
                this->index=0;
                this-> ztmp=0;
	}

	
	void FDPoissonBoltzmann_ICCG_PBC::initpciccg(std::vector<double> &ldiag,
                std::vector<double> &epsCgrid,
		std::vector<double> &epsIgrid,
                std::vector<double> &epsJgrid,
                std::vector<double> &epsKgrid) {
		
		
		index = 0;
               // double term1,term2,term3,term4,term5,term6;
                
		//construct ldiag once 
		for (int k=0; k < GPZ; ++k) {
			for (int j=0; j < GPY; ++j) {
				for (int i=0; i < GPX; ++i) { 
					index = (k) * GPXGPY + (j) * GPX + i;
					ldiag[index] = epsCgrid[index];

					if (i != 0) ldiag[index]     = ldiag[index] - pow(epsIgrid[index-1],2)/ldiag[index-1];
					if (i == GPX-1) ldiag[index] = ldiag[index] - pow(epsIgrid[index],2)/ldiag[index-GPX+1];					
					if (j != 0) ldiag[index]     = ldiag[index] - pow(epsJgrid[index-GPX],2)/ldiag[index-GPX];
					if (j == GPY-1) ldiag[index] = ldiag[index] - pow(epsJgrid[index],2)/ldiag[index-GPXGPY+GPX];
					if (k != 0) ldiag[index]     = ldiag[index] - pow(epsKgrid[index-GPXGPY],2)/ldiag[index-GPXGPY];
					if (k == GPZ-1) ldiag[index] = ldiag[index] - pow(epsKgrid[index],2)/ldiag[index-GPXGPYGPZ+GPXGPY];

                             /*           if (fabs(epsIgrid[index-1]) < tinynum   ){
                                        std::cout << "# !!! epsIgrid[index-1] "  << epsIgrid[index-1] << endl;
                                        term1=0.0;
                                        }
                                        else{
                                        term1=pow(epsIgrid[index-1],2)/ldiag[index-1];
                                        }
                                        if (fabs(epsIgrid[index]) < tinynum   ){
                                        std::cout << "# !!! epsIgrid[index] "  << epsIgrid[index]<< endl;
                                        term2=0.0;
                                        }
                                        else{
                                        term2=pow(epsIgrid[index],2)/ldiag[index-GPX+1];
                                        }
                                        if (fabs(epsJgrid[index-GPX]) < tinynum   ){
                                        std::cout << "# !!! epsJgrid[index-GPX] "  << epsJgrid[index-GPX]<< endl;
                                        term3=0.0;
                                        }
                                        else{
                                        term3=pow(epsJgrid[index-GPX],2)/ldiag[index-GPX];
                                        }
                                        if (fabs(epsJgrid[index]) < tinynum   ){
                                        std::cout << "# !!! epsJgrid[index] "  << epsJgrid[index]<< endl;
                                        term4=0.0;
                                        }
                                        else{
                                        term4=pow(epsJgrid[index],2)/ldiag[index-GPXGPY+GPX];
                                        }
                                        if (fabs(epsKgrid[index-GPXGPY]) < tinynum   ){
                                        std::cout << "# !!! epsKgrid[index-GPXGPY] "  << epsKgrid[index-GPXGPY]<< endl;
                                        term5=0.0;
                                        }
                                        else{
                                        term5=pow(epsKgrid[index-GPXGPY],2)/ldiag[index-GPXGPY];
                                        }
                                        if (fabs(epsKgrid[index]) < tinynum   ){
                                        std::cout << "# !!! epsKgrid[index] "  << epsKgrid[index]<< endl;
                                        term6=0.0;
                                        }
                                        else{
                                        term6=pow(epsKgrid[index],2)/ldiag[index-GPXGPYGPZ+GPXGPY];
                                        } 


                                        if (i != 0) ldiag[index]     = ldiag[index] - term1;
					if (i == GPX-1) ldiag[index] = ldiag[index] - term2;
					if (j != 0) ldiag[index]     = ldiag[index] - term3;
					if (j == GPY-1) ldiag[index] = ldiag[index] - term4;
					if (k != 0) ldiag[index]     = ldiag[index] - term5;
					if (k == GPZ-1) ldiag[index] = ldiag[index] - term6;
*/


				}
			}
		} 
		

		
	}
	
	
	void FDPoissonBoltzmann_ICCG_PBC::pciccg(std::vector<double> &zvec, std::vector<double>&ldiag,
                std::vector<double> &rhogrid,
		std::vector<double> &epsIgrid,
                std::vector<double> &epsJgrid,
                std::vector<double> &epsKgrid) {
				

		
	//int sizeomp = rhogrid.length;
	  int sizeomp = GPXGPYGPZ;

		//***** FORWARD SUBSTITUTION: INVERSE(L)*RHO - PERIODIC GRID ALONG X,Y,Z
		int k=0, j=0, i=0, index=0;
		double ztmp = 0.0;

		{		 
		 for (k=0; k < GPZ; ++k) {		 
			for (j=0; j < GPY; ++j) {				
				for (i=0; i < GPX; ++i) { 
					index = (k) * GPXGPY + (j) * GPX + i;
					ztmp = rhogrid[index];					
					if (i != 0)     ztmp += epsIgrid[index - 1] * zvec[index-1];
					if (i == GPX-1) ztmp += epsIgrid[index] * zvec[index-GPX+1];
					if (j != 0)     ztmp += epsJgrid[index-GPX] * zvec[index-GPX];
					if (j == GPY-1) ztmp += epsJgrid[index] * zvec[index-GPXGPY+GPX];
					if (k != 0)     ztmp += epsKgrid[index-GPXGPY] * zvec[index-GPXGPY];
					if (k == GPZ-1) ztmp += epsKgrid[index] * zvec[index-GPXGPYGPZ+GPXGPY];					
					 zvec[index] = ztmp/ldiag[index];				    
					
				}
			}
		 }
		}
		

		
		for (i=0; i < sizeomp-1; ++i) zvec[i] *= ldiag[i];

	
		
		//C***** BACK SUBSITUTION: INVERSE(TRANSPOSE(L))*ZVEC - PERIODIC GRID ALONG X,Y,
		ztmp = 0.0;
		
		{		 
		 for (k=GPZ-1; k >=0; --k) {
			for (j=GPY-1; j >=0; --j) {
				for (i=GPX-1; i >=0; --i) {
					index = (k) * GPXGPY + (j) * GPX + i;
					ztmp = zvec[index];
					
					if (i != GPX-1) ztmp += epsIgrid[index] * zvec[index+1];
					if (i == 0) ztmp += epsIgrid[index] * zvec[index+GPX-1];
					if (j != GPY-1) ztmp += epsJgrid[index] * zvec[index+GPX];
					if (j == 0) ztmp += epsJgrid[index] * zvec[index+GPXGPY-GPX];
					if (k != GPZ-1) ztmp += epsKgrid[index] * zvec[index+GPXGPY];
					if (k == 0) ztmp += epsKgrid[index] * zvec[index+GPXGPYGPZ-GPXGPY];
					zvec[index] = ztmp/ldiag[index]; 
				
				}
			}
		 }
		}
			 
		
	} //end pciccg

	

	
	 void FDPoissonBoltzmann_ICCG_PBC::gqact(std::vector<double> &zvec, std::vector<double> &pvec,
	         std::vector<double> &epsCgrid,
                 std::vector<double> &epsIgrid,
                 std::vector<double> &epsJgrid,
                 std::vector<double> &epsKgrid) {
		
		
		
		int k=0, j=0, i=0, index=0;
		double ztmp = 0.0;
	
		{		
		 for (k=0; k < GPZ; ++k) {
			for (j=0; j < GPY; ++j) {
				for (i=0; i < GPX; ++i) { 
					
					index = (k) * GPXGPY + (j) * GPX + i;
					ztmp  = epsCgrid[index] * pvec[index];

					if (i != GPX-1) ztmp -= epsIgrid[index] * pvec[index+1];
					else            ztmp -= epsIgrid[index] * pvec[index-GPX+1];
					
					if (i != 0)     ztmp -= epsIgrid[index-1] * pvec[index-1];
					else            ztmp -= epsIgrid[index+GPX-1] * pvec[index+GPX-1];
					
					if (j != GPY-1) ztmp -= epsJgrid[index] * pvec[index+GPX];
					else            ztmp -= epsJgrid[index] * pvec[index-GPXGPY+GPX];
					
					if (j != 0)     ztmp -= epsJgrid[index-GPX] * pvec[index-GPX];
					else            ztmp -= epsJgrid[index+GPXGPY-GPX] * pvec[index+GPXGPY-GPX];
					
					if (k != GPZ-1) ztmp -= epsKgrid[index] * pvec[index+GPXGPY];
					else            ztmp -= epsKgrid[index] * pvec[index-GPXGPYGPZ+GPXGPY];
					
					if (k != 0)     ztmp -= epsKgrid[index-GPXGPY] * pvec[index-GPXGPY];
					else            ztmp -= epsKgrid[index+GPXGPYGPZ-GPXGPY] * pvec[index+GPXGPYGPZ-GPXGPY];
					zvec[index] = ztmp;
				}
			}
		 }
		}
		
		
  
		
		
	}
	
	
