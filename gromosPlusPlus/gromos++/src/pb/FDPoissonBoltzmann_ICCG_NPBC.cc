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

// pb_FDPoissonBoltzmann_ICCG_NPBC.cc
#include "FDPoissonBoltzmann_ICCG_NPBC.h"

#include <cmath>
#include <cstdlib>
#include <cassert>
#include <vector>

using pb::FDPoissonBoltzmann_ICCG_NPBC;

FDPoissonBoltzmann_ICCG_NPBC::FDPoissonBoltzmann_ICCG_NPBC(int GPX, int GPY, int GPZ){

	
	
		this->GPX = GPX;
		this->GPY = GPY;
		this->GPZ = GPZ;
		this->GPXGPY = GPX*GPY;
		this->GPXGPYGPZ = GPX*GPY*GPZ;
                this->index=0;
          
	}

	 void FDPoissonBoltzmann_ICCG_NPBC::initpciccg(std::vector<double> &ldiag,
                 std::vector<double> & epsCgrid,
	         std::vector<double> & epsIgrid,
                 std::vector<double> & epsJgrid,
                 std::vector<double> & epsKgrid) {


        //     std::cout << "# GPX " << GPX << " GPXGPY " << GPXGPY << " GPXGPYPGPZ " << GPXGPYGPZ  << endl;


		double alpha = 0.0;
		//make the preconditioner...
		ldiag[0] = epsCgrid[0];
		for (int i = 1; i < GPX; ++i) {
               //     std::cout << "#----------" << endl;
                //    std::cout << "# loop over 1..GPX, fill ldiag[i], i = " << i << endl;
		 ldiag[i] = epsCgrid[i] 
		            - (pow(epsIgrid[i-1],2)
                 + (alpha * epsJgrid[i-1] * epsIgrid[i-1])
		 + (alpha * epsKgrid[i-1] * epsIgrid[i-1]))/ldiag[i-1];
                // std::cout << "# loop over 1..GPX, fill ldiag[i], ldiag[i] = " << ldiag[i] << endl;
		}

		for (int i = GPX; i < GPXGPY; ++i) {
                     // std::cout << "#----------" << endl;
                     // std::cout << "# loop over GPX..GPXGPY, fill ldiag[i], i = " << i << endl;
                     // std::cout << "# division by ldiag[i-1] which is = " <<  ldiag[i-1] << endl;
                     // std::cout << "# division by ldiag[i-GPX] which is = " <<  ldiag[i-1] << endl;
		 ldiag[i] = epsCgrid[i] 
		            - pow(epsIgrid[i-1],2)/ldiag[i-1] 
		            - (pow(epsJgrid[i-GPX],2)
		            + (alpha * epsIgrid[i-GPX] * epsJgrid[i-GPX])
		            + (alpha * epsKgrid[i-GPX] * epsJgrid[i-GPX]))/ldiag[i-GPX];
               //  std::cout << "# loop over GPX..GPXGPY, fill ldiag[i], ldiag[i] = " << ldiag[i] << endl;
		}

		for(int i = GPXGPY; i < GPXGPYGPZ; ++i) {
                 // std::cout << "#----------" << endl;
                 // std::cout << "# loop over GPXGPY..GPXGPYGPZ, fill ldiag[i], i = " << i << endl;
                 // std::cout << "# division by ldiag[i-1] which is = " <<  ldiag[i-1] << endl;
                //  std::cout << "# division by ldiag[i-GPX] which is = " <<  ldiag[i-GPX] << endl;
                //  std::cout << "# division by ldiag[i-GPXGY] which is = " <<  ldiag[i-GPXGPY] << endl;


                //  double term1=epsCgrid[i];
                 // double term2=0.0;
                 // double term3=0.0;
                //  double term4=0.0;

                  ldiag[i] = epsCgrid[i]
		            - pow(epsIgrid[i-1],2)/ldiag[i-1] 
		            -  pow(epsJgrid[i-GPX],2)/ldiag[i-GPX]
		            - (pow(epsKgrid[i-GPXGPY],2)
		            + (alpha * epsIgrid[i-GPXGPY] * epsKgrid[i-GPXGPY])
		            + (alpha * epsJgrid[i-GPXGPY] * epsKgrid[i-GPXGPY]))/ldiag[i-GPXGPY];
                 
/*
                      std::cout << "# loop over GPXGPY..GPXGPYGPZ, fill ldiag[i], ldiag[i] = " << ldiag[i] << endl;
               
                 
                 if ( fabs(epsIgrid[i-1]) < tinynum ){
                     std::cout << "# !!! epsIgrid[i] " << epsIgrid[i] << endl;
                     term2=0.0;
                 }
                 else{
                     term2= - pow(epsIgrid[i-1],2)/ldiag[i-1];
                 }
                 if ( fabs(epsJgrid[i-GPX]) < tinynum ){
                     std::cout << "# !!! epsJgrid[i-GPX] " << epsJgrid[i-GPX] << endl;
                     term3=0.0;
                 }
                 else{
                     term3=  -  pow(epsJgrid[i-GPX],2)/ldiag[i-GPX];
                 }
                 if ( fabs(epsKgrid[i-GPXGPY]) < tinynum ){
                     std::cout << "# !!! epsKgrid[i-GPXGPY] " << epsKgrid[i-GPXGPY] << endl;
                     term4=0.0;
                 }
                 else{
                     term4=  - (pow(epsKgrid[i-GPXGPY],2)
		            + (alpha * epsIgrid[i-GPXGPY] * epsKgrid[i-GPXGPY])
		            + (alpha * epsJgrid[i-GPXGPY] * epsKgrid[i-GPXGPY]))/ldiag[i-GPXGPY];
                 }


                 ldiag[i]=term1+term2+term3+term4;
 */
                     }
 
		

	}
	
	
        void FDPoissonBoltzmann_ICCG_NPBC::pciccg(std::vector<double> & zvec, std::vector<double> & ldiag,
                std::vector<double> & rhogrid,
	        std::vector<double> & epsIgrid,
                std::vector<double> & epsJgrid,
                std::vector<double> & epsKgrid) {

		//***** FORWARD SUBSTITUTION: INVERSE(L)*RHO - NO PERIODICITY
		
		//handle implicit boundary point first
		
		zvec[0] = rhogrid[0] / ldiag[0];
		
		for (int i = 1; i < GPX; ++i) {          
			zvec[i] = ( rhogrid[i] + epsIgrid[i-1]*zvec[i-1] ) / ldiag[i];       
		}
	
		for (int i = GPX; i < GPXGPY; ++i) {        
			zvec[i] = ( rhogrid[i] + (epsIgrid[i-1]*zvec[i-1]) + (epsJgrid[i-GPX]*zvec[i-GPX])) / ldiag[i];        
		} 
		
		for(int i = GPXGPY; i < GPXGPYGPZ; ++i) {
			zvec[i] = ( rhogrid[i] + epsIgrid[i-1]*zvec[i-1] + epsJgrid[i-GPX]*zvec[i-GPX] + epsKgrid[i-GPXGPY]*zvec[i-GPXGPY] ) / ldiag[i];        
		}
		
		
		
		//***** Z = Z * L; add some OpenMP directives
		
		for (int i=0; i < GPXGPYGPZ-1; ++i) zvec[i] *= ldiag[i];

		//***** BACK SUBSITUTION: INVERSE(TRANSPOSE(L))*ZVEC - NO PERIODICITY
		
		
		for (int i = GPXGPYGPZ-2; i >= GPXGPYGPZ-GPX; --i) {
			zvec[i] = (zvec[i] + epsIgrid[i]*zvec[i+1]) / ldiag[i];      
		}
		
		for (int i = GPXGPYGPZ-GPX-1; i >= GPXGPYGPZ-GPXGPY; --i) {
			zvec[i] = (zvec[i] + epsIgrid[i]*zvec[i+1]+ epsJgrid[i]*zvec[i+GPX]) / ldiag[i];       
		}
	
		for (int i = GPXGPYGPZ-GPXGPY-1; i >= 0; --i) { 
			zvec[i] = (zvec[i] + epsIgrid[i]*zvec[i+1] + epsJgrid[i]*zvec[i+GPX] + epsKgrid[i]*zvec[i+GPXGPY]) / ldiag[i];      
		}
        
		
		
	} //end pciccg
	
	
       void FDPoissonBoltzmann_ICCG_NPBC::gqact(std::vector<double> & zvec,std::vector<double> & pvec,
					 std::vector<double> & epsCgrid,
                                         std::vector<double> & epsIgrid,
                                         std::vector<double> & epsJgrid,
                                         std::vector<double> & epsKgrid) {

		//long time = System.currentTimeMillis();
		
		//diagonal
				for (int i=0; i < GPXGPYGPZ; ++i) zvec[i] = epsCgrid[i] * pvec[i];
		//  First super diagonal
				for (int i=0; i < GPXGPYGPZ-1; ++i) zvec[i] = zvec[i] - epsIgrid[i] * pvec[i+1];
		//  Second super diagonal
			for (int i=0; i < GPXGPYGPZ-GPX; ++i) zvec[i] = zvec[i] - epsJgrid[i] * pvec[i+GPX];
		//  Third super diagonal
				for (int i=0; i < GPXGPYGPZ-GPXGPY; ++i) zvec[i] = zvec[i] - epsKgrid[i] * pvec[i+GPXGPY];
		//  First sub diagonal
			for (int i=1; i < GPXGPYGPZ; ++i) zvec[i] = zvec[i] - epsIgrid[i-1] * pvec[i-1];
		//  Second sub diagonal
			for (int i=GPX; i < GPXGPYGPZ; ++i) zvec[i] = zvec[i] - epsJgrid[i-GPX] * pvec[i-GPX];
		//  Third sub diagonal
				for (int i=GPXGPY; i < GPXGPYGPZ; ++i) zvec[i] = zvec[i] - epsKgrid[i-GPXGPY] * pvec[i-GPXGPY];
        
        
	}
	
	
