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

// pb_Ewald_edir.cc
#include "Ewald_edir.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <complex>
#include <vector>

#include "PB_Parameters.h"
#include "../fit/PositionUtils.h"
#include "../utils/AtomSpecifier.h"
#include "../gcore/System.h"
#include "../gcore/Box.h"
#include "../bound/Boundary.h"
#include "../gromos/Exception.h"

using pb::Ewald_edir;
using namespace std;

 Ewald_edir::Ewald_edir(bound::Boundary & pbc, utils::AtomSpecifier atoms,
			double realcut, double tolerance, int KX_in, int KY_in, int KZ_in, ofstream &os)
 //bool exinterx)
 : ppp(os){
     
    // this-> excluded_interx=exinterx;

     this-> pbc = &pbc;
     this-> atoms=atoms;
   
     this-> realcut=realcut;
     this-> tolerance=tolerance;
     
     this-> kx=KX_in;
     this-> ky=KY_in;
     this-> kz=KZ_in;

     int nx = kx+1;
     int ny = ky+1;
     int nz = kz+1;


    // the maximum: kmax
     this->kmax=nx;
     if (ny>kmax){kmax=ny;}
     if (nz>kmax){kmax=nz;}


     
     
   // resize the eir: 0: kmax; 1: atoms.size(); 2: 3
     
     
     eir.resize(kmax);
     for (int i=0; i<kmax; i++){
         eir[i].resize(atoms.size());
            for (unsigned int j=0; j<atoms.size(); j++){
                 eir[i][j].resize(3);
             }}
     
  /*   eir_slow.resize(kmax*2+1);
     for (int i=-(kmax-1); i<kmax; i++){
         eir_slow[i].resize(atoms.size());
            for (int j=0; j<atoms.size(); j++){
                 eir_slow[i][j].resize(3);
             }}*/





     this-> thebox = pbc.sys().box();

    

     // for convenience, a box array:
     this->box[0]=thebox.K().abs();
     this->box[1]=thebox.L().abs();
     this->box[2]=thebox.M().abs();



     std::cout<<"# box x, y, z: " << box[0] << " " << box[1] << " " << box[2]  << endl;

     // exit if the system is not rectangular

     try{
         if ( pbc.type() != 'r'){
                throw gromos::Exception("Ewald_edir","System does not have rectangular spatial boundary conditions ...");
              }// end of if
	 } // end of try
                                catch (const gromos::Exception &e){
                                         cerr << e.what() << endl;
                                         exit(1);
                                }


     
 }

   complex<double> Ewald_edir::scalecomplexnum(   double r,
				                       complex<double> c) {

		

			 complex<double>           ret ( r * real(c),
			                              r * imag(c));

			return ret;
		}

  complex<double> Ewald_edir::multiplycomplexnums(
                                 complex<double> a,
                                 complex<double> b) {
                      

                         complex<double> ret (real(a) * real(b) - imag(a) * imag(b),
                                              real(a) * imag(b) + imag(a) * real(b));

                        return ret;
                }




  complex<double> Ewald_edir::multiplycomplexnums_firstconj(
                                 complex<double> a,
                                 complex<double> b) {


                         complex<double> ret (real(a) * real(b) + imag(a) * imag(b),
                                              real(a) * imag(b) - imag(a) * real(b));

                        return ret;
                }


  complex<double> Ewald_edir::multiplycomplexnums_bothconj(
                                 complex<double> a,
                                 complex<double> b) {


                         complex<double> ret (real(a) * real(b) - imag(a) * imag(b),
                                              -real(a) * imag(b) - imag(a) * real(b));

                        return ret;
                }


   complex<double> Ewald_edir::multiplycomplexconjugate(
				                        complex<double> a,
				                        complex<double> b) {

			

			 complex<double>  ret  ( real(a) * real(b) + imag(a) * imag(b),
			                        imag(a) * real(b) - real(a) * imag(b));

			return ret;
		}

	


    void Ewald_edir::calcenergy(double (& energy)[15]) {



		// Ewald coefficient beta
		double ewaldcoeff = calc_ewaldcoeff();
                cout << "# ewaldcoeff " << ewaldcoeff << endl;
                
		// real-space energy
		energy[0] = rspaceEwald(ewaldcoeff);
		// reciprocal (k-)space Ewald energy
		energy[1] = kspaceEwald(ewaldcoeff);
		// self-term
		energy[2] = -minusA3timesStildeSquare(ewaldcoeff);
		// charge correction
		energy[3] = -minusA1timesSSquare(ewaldcoeff);
		// sum of energy components
		energy[4] = energy[0] + energy[1] + energy[2] + energy[3];

                
                // self-term with sumq * xiew/l
                energy[5] = self_other();
                // coulombic interx between non-excluded atoms
                energy[6] = coulomb_non_excluded_atoms();
                energy[7] = XIEWcorr();
               
                energy[8] =  A2timesStildeSquare(ewaldcoeff);
                energy[9] = A1timesStildeSquare(ewaldcoeff);
                
                // the EA of Phil : EA = (A1 S2 - A1 S~2 - A2 S~2 )
                energy[10] = energy[3] - energy[9] - energy[8] ;
                // the Eslf of Phil, now via sum of A1,A2,A3
                energy[11] = energy[9] + energy[8] + energy[2] ;

                // sum with Eslf, 11
                energy[12] = energy[0] + energy[1] + energy[10] + energy[11];
                // sum with Eslf, 5
                energy[13] = energy[0] + energy[1] + energy[10] + energy[5];


                energy[14]=0.0;
               // energy[14]=kspaceEwald_slow(ewaldcoeff);

	}


    double Ewald_edir::rspaceEwald(
			double ewaldcoeff) {

		double energy = 0.0;

		for (unsigned int i=0; i < atoms.size()-1; ++i) {
			for (unsigned int j=i+1; j < atoms.size(); ++j) {
                      
                            

                        // get the distance between the atoms
                          gmath::Vec jimage = (*pbc).nearestImage(atoms.pos(i),atoms.pos(j), thebox);
                          gmath::Vec dd = atoms.pos(i)-jimage;
                          double distance=dd.abs();


                          //bool excluded_atoms;


                          //int im=atoms.mol(i);
                          //int ia=atoms.atom(i);
                          //int jm=atoms.mol(j);
                          //int ja=atoms.atom(j);

                      


                          //gcore::System *sss= atoms.sys();
                          
                         //excluded_atoms = false;

                        //  if (im==jm && sss->mol(im).topology().atom(ia).exclusion().contains(ja)){
                        //      excluded_atoms=true;
                       //   }
                          
                         
                          
                        //  if (!excluded_atoms){
                  	     // energy += ((atoms.charge(i) * atoms.charge(j))/distance) * erfcapp(distance*ewaldcoeff);
                              energy += ((atoms.charge(i) * atoms.charge(j))/distance) *  (erfc(distance*ewaldcoeff)  ); //+ 2.837297/box[0]);
                             // cout << "# atoms i,j: " <<  i << " " << j << endl;
                             // cout << "# charges of atoms i,j: " << atoms.charge(i) << " " << atoms.charge(j) << endl;
                             // cout << "# distance " << distance << endl;
                            //  cout << "# energy now " << energy << endl;
                         // }
                         // else{

                          //    if (excluded_interx){
                               // they may interact, but remove the Coulomb contribution
                            //  energy += ((atoms.charge(i) * atoms.charge(j))/distance) * (erfcapp(distance*ewaldcoeff) - 1.0);
                         //    energy += ((atoms.charge(i) * atoms.charge(j))/distance) * erfc(distance*ewaldcoeff);
                         //     }
                          //}

                        }
		}

		return (energy *  ppp.getFPEPSI() );
	}



	//takes real-space cutoff and tolerance as input
	/* Real space tolerance for Ewald, determines   */
	/* the real/reciprocal space relative weight    */
	double Ewald_edir::calc_ewaldcoeff() {
		double x=5,low,high;
		int n,i=0;

		do {
			i++;
			x*=2;
		} while (erfc(x*realcut) > tolerance);

		n=i+60; /* search tolerance is 2^-60 */
		low=0;
		high=x;
		for(i=0;i<n;i++) {
			x=(low+high)/2.0;
			if (erfc(x*realcut) > tolerance) {//(jpb_solver_Math.erfc(x*rc) > dtol) {
				low=x;
			}
			else {
				high=x;
			}
		}

                std::cout << "calc_ewaldcoef alpha = a^(-1) = " << x << " [for: gamma(a^(-1)  * r) = pi^(-3/2) exp[- (a^(-1) r)^2]]\n";
		return x;
	}

	void Ewald_edir::tabulate_eir(double (& lll)[3]) {

            
	unsigned int  i;
        int j,m;

              try{
		if (kmax < 1) {
			       throw gromos::Exception("tabulate_eir","kmax < 1 ...");
                               } 
               
                else{

		//make structure factor for each atom...

	

			for(i=0; i<atoms.size(); i++) {
				for(m=0; m<3; m++) {
					eir[0][i][m] = complex<double>(1.0,0.0);
					
				}
				for(m=0; m<3; m++) {
					eir[1][i][m] = complex<double> (cos(atoms.pos(i)[m]*lll[m]),
                                                                         sin(atoms.pos(i)[m]*lll[m]));
			
				}
				for(j=2; j<kmax; j++) {
					for(m=0; m<3; m++) {
						eir[j][i][m] = multiplycomplexnums(eir[j-1][i][m],eir[1][i][m]);
					}
				}
			}

                }// end of if
		  } // end of try
                                catch (const gromos::Exception &e){
                                         cerr << e.what() << endl;
                                         exit(1);
                                }
	}

	double  Ewald_edir::kspaceEwald(
			double ewaldcoeff) {


                     int nx = kx+1;
                     int ny = ky+1;
                     int nz = kz+1;

		
		std::vector< complex <double> > tab_xy;
		std::vector< complex <double> > tab_qxyz;

                
		double factor=-1.0/(4*ewaldcoeff*ewaldcoeff);

                double lll[3];

	
			for (int i=0; i < kmax; ++i) {
				for (unsigned int ii=0; ii < atoms.size(); ++ii) {
					for (int iii=0; iii < 3; ++iii) {
						 eir[i][ii][iii]= complex<double> (0.0, 0.0);
                                              					}
				}
			}
		

		//tabulation ComplexDouble arrays...
		tab_xy.resize(atoms.size());
		tab_qxyz.resize(atoms.size());

		calc_lll(lll);

            //    cout << "# lll " << lll[0] << " " << lll[1] << " " << lll[2] << endl;


		/* make tables for the structure factor parts */
		tabulate_eir(lll);


		int  ix,iy,iz;
                unsigned int n;
		double tmp,cs,ss,ak,akv,mx,my,mz,m2;
		double energy=0.0;
		int lowiy=0;
		int lowiz=1;

                double a2sumtest = 0;
		for(ix=0;ix<nx;ix++) {
			mx=ix*lll[0];
			for(iy=lowiy;iy<ny;iy++) {
				my=iy*lll[1];
				if (iy>=0) {
				
					for(n=0;n<atoms.size();n++) {
						tab_xy[n]=multiplycomplexnums(eir[ix][n][0],eir[iy][n][1]);
					}
				}
				else {
				
					for(n=0;n<atoms.size();n++) {
						tab_xy[n]=multiplycomplexconjugate(eir[ix][n][0],eir[-iy][n][1]);
					}
				}
				for(iz=lowiz;iz<nz;iz++) {
					mz=iz*lll[2];
					m2=mx*mx+my*my+mz*mz;
					ak=exp(m2*factor)/m2;
					akv=2.0*ak*(1.0/m2-factor);
					if(iz>=0) {
					
						for(n=0;n<atoms.size();n++) {
							tab_qxyz[n]=scalecomplexnum(atoms.charge(n),
                                                                    multiplycomplexnums(tab_xy[n],eir[iz][n][2]));
						}
					}
					else {
				
						for(n=0;n<atoms.size();n++) {
							tab_qxyz[n]=scalecomplexnum(atoms.charge(n),
                                                                    multiplycomplexconjugate(tab_xy[n],eir[-iz][n][2]));
						}
					}
					cs=ss=0;
					
					for(n=0;n<atoms.size();n++) {
						cs+=tab_qxyz[n].real();
						ss+=tab_qxyz[n].imag();
					}
					energy+=ak*(cs*cs+ss*ss);
					tmp=akv*(cs*cs+ss*ss);

                                        a2sumtest +=ak;
					lowiz=1-nz;
				}
				lowiy=1-ny;
			}
		}
               
		tmp=4.0*ppp.getPI()/(box[0]*box[1]*box[2]) * ppp.getFPEPSI();

                a2sumtest*=tmp;
		energy*=tmp;
                std::cout << "a2sumtest " << a2sumtest << "\n";
	

		return energy;
	}

	double  Ewald_edir::kspaceEwald_slow(
			double ewaldcoeff) {


                     int nx = kx+1;
                     int ny = ky+1;
                     int nz = kz+1;


		std::vector< complex <double> > tab_xy;
		std::vector< complex <double> > tab_qxyz;


		double factor=-1.0/(4*ewaldcoeff*ewaldcoeff);

                double lll[3];


			for (int i=0; i < kmax; ++i) {
				for (unsigned int ii=0; ii < atoms.size(); ++ii) {
					for (int iii=0; iii < 3; ++iii) {
						 eir[i][ii][iii]= complex<double> (0.0, 0.0);
                                              					}
				}
			}


		//tabulation ComplexDouble arrays...
		tab_xy.resize(atoms.size());
		tab_qxyz.resize(atoms.size());

		calc_lll(lll);

            //    cout << "# lll " << lll[0] << " " << lll[1] << " " << lll[2] << endl;


		/* make tables for the structure factor parts */
		tabulate_eir(lll);


		int  ix,iy,iz;
                unsigned int n;
		double tmp,cs,ss,ak,akv,mx,my,mz,m2;
		double energy=0.0;

                // int lowix=0;
                //int lowiy=0;
		//int lowiz=1;

                int lowix=1-nx;
                int lowiy=1-ny;
		int lowiz=1-nz;

                
                double a2sumtest = 0;
		for(ix=lowix;ix<nx;ix++) {
			mx=ix*lll[0];
			for(iy=lowiy;iy<ny;iy++) {
				my=iy*lll[1];


                                if (ix>=0){

				if (iy>=0) {

					for(n=0;n<atoms.size();n++) {
						tab_xy[n]=multiplycomplexnums(eir[ix][n][0],eir[iy][n][1]);
					}
				}
				else {

					for(n=0;n<atoms.size();n++) {
						tab_xy[n]=multiplycomplexconjugate(eir[ix][n][0],eir[-iy][n][1]);
					}
				}
                                }
                                else{
                                  if (iy>=0) {

					for(n=0;n<atoms.size();n++) {
						tab_xy[n]=multiplycomplexnums_firstconj(eir[-ix][n][0],eir[iy][n][1]);
					}
				}
				else {

					for(n=0;n<atoms.size();n++) {
						tab_xy[n]=multiplycomplexnums_bothconj(eir[-ix][n][0],eir[-iy][n][1]);
					}
				}
                                }

				for(iz=lowiz;iz<nz;iz++) {


                                    if ( ! (ix==0 && iy==0 && iz==0) ) {


					mz=iz*lll[2];
					m2=mx*mx+my*my+mz*mz;
					ak=exp(m2*factor)/m2;
					akv=2.0*ak*(1.0/m2-factor);
					if(iz>=0) {

						for(n=0;n<atoms.size();n++) {
							tab_qxyz[n]=scalecomplexnum(atoms.charge(n),
                                                                    multiplycomplexnums(tab_xy[n],eir[iz][n][2]));
						}
					}
					else {

						for(n=0;n<atoms.size();n++) {
							tab_qxyz[n]=scalecomplexnum(atoms.charge(n),
                                                                    multiplycomplexconjugate(tab_xy[n],eir[-iz][n][2]));
						}
					}
					cs=ss=0;

					for(n=0;n<atoms.size();n++) {
						cs+=tab_qxyz[n].real();
						ss+=tab_qxyz[n].imag();
					}
					energy+=ak*(cs*cs+ss*ss);
					tmp=akv*(cs*cs+ss*ss);

                                        a2sumtest +=ak;

                                    }


				//	lowiz=1-nz;
				}
				//lowiy=1-ny;
			}
		}

		tmp=4.0*ppp.getPI()/(box[0]*box[1]*box[2]) * ppp.getFPEPSI() / 2.0 ;

                a2sumtest*=tmp;
		energy*=tmp;
                std::cout << "a2sumtest " << a2sumtest << "\n";


		return energy;
	}


        double  Ewald_edir::A2timesStildeSquare(
			double ewaldcoeff) {


                     int nx = kx;
                     int ny = ky;
                     int nz = kz;


		//std::vector< complex <double> > tab_xy;
		//std::vector< complex <double> > tab_qxyz;

                double a2sum = 0.0;
		double factor=-1.0/(4*ewaldcoeff*ewaldcoeff);

                double lll[3];


		//	for (int i=0; i < kmax; ++i) {
		//		for (int ii=0; ii < atoms.size(); ++ii) {
		//			for (int iii=0; iii < 3; ++iii) {
		//				 eir[i][ii][iii]= complex<double> (0.0, 0.0);
                  //                            					}
		//		}
		//	}


		//tabulation ComplexDouble arrays...
		//tab_xy.resize(atoms.size());
		//tab_qxyz.resize(atoms.size());

		calc_lll(lll);

            //    cout << "# lll " << lll[0] << " " << lll[1] << " " << lll[2] << endl;


		/* make tables for the structure factor parts */
		//tabulate_eir(lll);


		int  ix,iy,iz;
		double tmp,ak,mx,my,mz,m2;
	



		for(ix=-nx;ix<=nx;ix++) {
			mx=ix*lll[0];
			for(iy=-ny;iy<=ny;iy++) {
				my=iy*lll[1];
				
				for(iz=-nz;iz<=nz;iz++) {
                                    
                                      if ( ! (ix==0 && iy==0 && iz==0) ){
                                    
					mz=iz*lll[2];
					m2=mx*mx+my*my+mz*mz;
					ak=exp(m2*factor)/m2;



                                      
					a2sum+=ak;
                                        }


				}
			
			}
		}

		tmp=4.0*ppp.getPI()/(box[0]*box[1]*box[2]) * ppp.getFPEPSI() / 2.0;


		a2sum*=tmp;


                // now get stildesquare
                double q2sum=0;
               for (unsigned int i=0; i < atoms.size();++i){
                    q2sum += atoms.charge(i) * atoms.charge(i);
               }
                std::cout << "a2sum " << a2sum << "\n";
		return (a2sum * q2sum);
	}








	double Ewald_edir::minusA3timesStildeSquare(
			double ewaldcoeff) {

		double q2_sum = 0.0;
		for (unsigned int i=0; i < atoms.size();++i){
                    q2_sum += atoms.charge(i) * atoms.charge(i);
                }
		return ewaldcoeff * ppp.getFPEPSI() * q2_sum/(sqrt(ppp.getPI()));

	}

	
	double Ewald_edir::minusA1timesSSquare(
			double ewaldcoeff) {
		double q_sum = 0.0;
		for (unsigned int i=0; i < atoms.size();++i) {q_sum += atoms.charge(i);};
	
		double vol = box[0]*box[1]*box[2];
		return((q_sum*q_sum * ppp.getPI() * ppp.getFPEPSI())/(2.0*vol*ewaldcoeff*ewaldcoeff));
	}

        double Ewald_edir::A1timesStildeSquare(
			double ewaldcoeff) {
		double q_sum = 0.0;
		for (unsigned int i=0; i < atoms.size();++i){ q_sum += (atoms.charge(i) * atoms.charge(i));};

		double vol = box[0]*box[1]*box[2];
		return( -1.0 * (q_sum * ppp.getPI() * ppp.getFPEPSI())/(2.0*vol*ewaldcoeff*ewaldcoeff));
	}




	void Ewald_edir::calc_lll(double (& lll) [3]) {
		for (int i=0; i < 3; ++i) lll[i] = 2.0 * ppp.getPI()/box[i];
	}


        double Ewald_edir::self_other(){
            double qsum = 0;
            for (unsigned int i=0; i < atoms.size();++i) qsum += atoms.charge(i)*atoms.charge(i);
            double res= ppp.getFPEPSI()/2.0 * qsum / box[0] * (ppp.get_xiew()) ;
            return res;
        }


      /*  double Ewald_edir::dip2(){
            // compute the squared dip moment w.r.t. the center of geom
            gmath::Vec cc;
            gmath::Vec dip;
            for (int i=0; i<3; ++i){
                dip[i]=0.0;
                cc[i]=0.0;
            }
         
            for (int i=0; i<atoms.size(); ++i){
            	for (int j=0; j<3; ++j){
               	 cc[j]=0.0;
           	 }
               for (int j=0; j<atoms.size(); ++j){
                   if ( j!=i){
                   cc+=atoms.pos(j);
                    }
                }
                 cc=cc/(atoms.size()-1.0);
                dip+=(atoms.pos(i)-cc)*atoms.charge(i);
              
            }
             std::cout << "@ dip " << dip[0] << " " << dip[1] << " " << dip[2] << "\n";
            double dipmag2 = dip.abs() * dip.abs();
            return dipmag2;
        }*/


   double Ewald_edir::coulomb_non_excluded_atoms(
		) {

		double energy = 0.0;

		for (unsigned int i=0; i < atoms.size()-1; ++i) {
			for (unsigned int j=i+1; j < atoms.size(); ++j) {



                        // get the distance between the atoms
                          gmath::Vec jimage = (*pbc).nearestImage(atoms.pos(i),atoms.pos(j), thebox);
                          gmath::Vec dd = atoms.pos(i)-jimage;
                          double distance=dd.abs();


                          //bool excluded_atoms;


                          //int im=atoms.mol(i);
                          //int ia=atoms.atom(i);
                          //int jm=atoms.mol(j);
                          //int ja=atoms.atom(j);




                          //gcore::System *sss= atoms.sys();

                          //excluded_atoms = false;

                        //  if (im==jm && sss->mol(im).topology().atom(ia).exclusion().contains(ja)){
                        //      excluded_atoms=true;
                        //  }

                          

                       //   if (!excluded_atoms){

                              energy += ((atoms.charge(i) * atoms.charge(j))/distance);
                           //   std::cout << "i " << i << " j " << j << " charges "  << atoms.charge(i)  << " " << atoms.charge(j)  << " dist " << distance << "\n";
     
                        //  }


                        }
		}

		return (energy *  ppp.getFPEPSI() );
	}




    double Ewald_edir::XIEWcorr(
		) {

		double energy = 0.0;

		for (unsigned int i=0; i < atoms.size()-1; ++i) {
			for (unsigned int j=i+1; j < atoms.size(); ++j) {



                        // get the distance between the atoms
                          gmath::Vec jimage = (*pbc).nearestImage(atoms.pos(i),atoms.pos(j), thebox);
                          gmath::Vec dd = atoms.pos(i)-jimage;
                          //double distance=dd.abs();


                          //bool excluded_atoms;


                          //int im=atoms.mol(i);
                          //int ia=atoms.atom(i);
                          //int jm=atoms.mol(j);
                          //int ja=atoms.atom(j);




                          //gcore::System *sss= atoms.sys();

                       //   excluded_atoms = false;

                        //  if (im==jm && sss->mol(im).topology().atom(ia).exclusion().contains(ja)){
                         //     excluded_atoms=true;
                         // }



                         // if (!excluded_atoms){

                              energy += ((atoms.charge(i) * atoms.charge(j))/box[0]) * (-1.0 * ppp.get_xiew());
                             // std::cout << "i " << i << " j " << j << " charges "  << atoms.charge(i)  << " " << atoms.charge(j)  << " dist " << distance << "\n";

                         // }


                        }
		}

		return (energy *  ppp.getFPEPSI() );
	}




 // actually we do not use this one
	double Ewald_edir::erfcapp(double X) {


		//    *******************************************************************
		//    ** APPROXIMATION TO THE COMPLEMENTARY ERROR FUNCTION             **
		//    **                                                               **
		//    ** REFERENCE:                                                    **
		//    **                                                               **
		//    ** ABRAMOWITZ AND STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS,    **
		//    **    NATIONAL BUREAU OF STANDARDS, FORMULA 7.1.26               **
		//    *******************************************************************

		double A1 = 0.254829592, A2 = -0.284496736;
		double A3 = 1.421413741, A4 = -1.453152027;
		double A5 = 1.061405429, P  =  0.3275911;

		double T  = 1.0 / ( 1.0 + P * X );
		double XSQ = X * X;

		double TP = T * ( A1 + T * ( A2 + T * ( A3 + T * ( A4 + T * A5 ) ) ) );

		double ERFC = TP * exp(-XSQ);

		return ERFC;
	}


