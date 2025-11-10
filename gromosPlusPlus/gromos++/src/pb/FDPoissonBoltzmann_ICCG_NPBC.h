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

// pb_FDPoissonBoltzmann_ICCG_NPBC.h

#ifndef INCLUDED_PB_FDPoissonBoltzmann_ICCG_NPBC
#define INCLUDED_PB_FDPoissonBoltzmann_ICCG_NPBC

#include <vector>
namespace pb{

class FDPoissonBoltzmann_ICCG_NPBC{

	
	
	int GPX;
	int GPY;
	int GPZ;
	int GPXGPY;
	int GPXGPYGPZ;
	
	int index;
	
public:
           //constructor
           FDPoissonBoltzmann_ICCG_NPBC(int GPX, int GPY, int GPZ);
	   // deconstructor
           ~FDPoissonBoltzmann_ICCG_NPBC(){}

           //methods
	
	
            void initpciccg(std::vector<double> &ldiag,
            std::vector<double> &epsCgrid,
            std::vector<double> &epsIgrid,
            std::vector<double> &epsJgrid,
            std::vector<double> &epsKgrid);

	
	
	
	    void pciccg(std::vector<double> &zvec, std::vector<double> &ldiag,
            std::vector<double> &rhogrid,
	    std::vector<double> &epsIgrid,
            std::vector<double> &epsJgrid,
            std::vector<double> &epsKgrid);

	
            void gqact(std::vector<double> &zvec, std::vector<double> &pvec,
            std::vector<double> &epsCgrid,
            std::vector<double> &epsIgrid,
            std::vector<double> &epsJgrid,
            std::vector<double> &epsKgrid);




}; //class
} //namespace
#endif
