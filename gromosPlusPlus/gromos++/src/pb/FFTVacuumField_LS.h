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

// pb_FFTVacuumField_LS.h

#ifndef INCLUDED_PB_FFTVacuumField_LS
#define INCLUDED_PB_FFTVacuumField_LS

#include <vector>

#include "FFTVacuumField.h"
#include "FFTChargeShapingFunction.h"
#include "PB_Parameters.h"

namespace pb{



class FFTVacuumField_LS: virtual public FFTVacuumField{
    
        FFTChargeShapingFunction *csfunc;

        int cstype;

	std::vector <int> ion_list1;
	std::vector <int> ion_list2;
	
	int ion_count1;
	int ion_count2;
	
	
	 /* arrays for storing configuration independent part of the vacuum field */

	std::vector <double> Vx1_store;
	std::vector <double> Vy1_store;
	std::vector <double> Vz1_store;
	std::vector <double> Vx2_store;
	std::vector <double> Vy2_store;
	std::vector <double> Vz2_store;
	
	int ngrdx;
	int ngrdy;
	int ngrdz;
	double dkx;
	double dky;
	double dkz;
	

public:
        //constructor
        FFTVacuumField_LS(utils::AtomSpecifier atoms, FFTGridType gt, FFTBoundaryCondition bc, ofstream &os);
       //deconstructor
        ~FFTVacuumField_LS(){}
       //methods

	
	
	void complexFromDouble(
			std::vector <double> & doubleArray, std::vector <double> & complexArray);
		
		
	 /* make lists of "big" and "small" atoms */
	 
	void makeatomlist(ofstream &os);
	
	
	void calcVacField(
			std::vector <double> &fldx_k,
                        std::vector <double> & fldy_k,
                        std::vector <double> & fldz_k,
			ofstream &os);

        void positionIndependentVacField(
			std::vector <double> &  destX, std::vector <double> &  destY, std::vector <double> &  destZ,
			double kax, double kay, double kaz,
			int nAlias, double alpha, ofstream &os);

        void recyclefield(
			  std::vector <double> &  destX, std::vector <double> & destY, std::vector <double> &  destZ, ofstream &os);

	void updateMultipliers(
			std::vector <double> &  multipliers,
			std::vector <int> & atomIndices, int numAtoms);

}; // class
} // namespace


#endif
		
