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

// pb_Ewald_edir.h


#ifndef INCLUDED_PB_Ewald_edir
#define INCLUDED_PB_Ewald_edir

#include <complex>

#include "PB_Parameters.h"
#include "../bound/Boundary.h"
#include "../gcore/Box.h"
#include "../utils/AtomSpecifier.h"


namespace pb{



class Ewald_edir{
    
bound::Boundary *pbc;
utils::AtomSpecifier atoms;
pb::PB_Parameters ppp;
gcore::Box thebox;

std::vector<std::vector<std::vector< complex<double>  > > >  eir;
//std::vector<std::vector<std::vector< complex<double>  > > >  eir_slow;
bool excluded_interx;

int kmax;
double realcut;
double tolerance;
double box[3];


int kx;
int ky;
int kz;


 public:
  // constructor

  Ewald_edir(bound::Boundary & pbc,utils::AtomSpecifier atoms,
             double realcut, double tolerance,  int KX_in, int KY_in, int KZ_in, ofstream &os);// bool exinterx);


   // deconstructor
  ~Ewald_edir(){}



  //methods


     complex<double> scalecomplexnum(       double r,
				            complex<double> c);

     complex<double> multiplycomplexnums(
                                 complex<double> a,
                                 complex<double> b);

     complex<double> multiplycomplexconjugate(
				complex<double> a,
				complex<double> b);


     complex<double> multiplycomplexnums_firstconj(
				complex<double> a,
				complex<double> b);

     complex<double> multiplycomplexnums_bothconj(
				complex<double> a,
				complex<double> b);
	
        void calcenergy(double (& energy)[15]);


        double rspaceEwald(double ewaldcoeff);
	
	
	double calc_ewaldcoeff();
	
	void tabulate_eir(double (& lll)[3]);
	
	double kspaceEwald(
			double ewaldcoeff);


	double kspaceEwald_slow(
			double ewaldcoeff);
        
	double A2timesStildeSquare(
			double ewaldcoeff);

	double minusA3timesStildeSquare(
			double ewaldcoeff);


	double A1timesStildeSquare(
			double ewaldcoeff);


        
        double self_other();
	 double coulomb_non_excluded_atoms();
    double  XIEWcorr();
    //double dip2();
	double minusA1timesSSquare(
			double ewaldcoeff);
	
	void calc_lll(double (& lll)[3]);


        // actually we do not use this one
	double erfcapp(double X);




}; // class
} // namespace


#endif
