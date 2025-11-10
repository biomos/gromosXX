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

// pb_PB_Parameters.cc
#include "PB_Parameters.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>

#include "../gmath/Physics.h"


using pb::PB_Parameters;
using namespace std;

PB_Parameters::PB_Parameters(double epssolvent, ofstream &os){
  os << "# initialising PB_Parameters using epsilon of " << epssolvent << " for adjusted boundary conditions" << endl;
  initParams(epssolvent, os);
}




PB_Parameters::PB_Parameters(ofstream &os){
  os << "# initialising PB_Parameters using default epsilon of 78.4 for adjusted boundary conditions" << endl;
  initParams(78.4, os);// DEFAULT ADJUSTED BOUNDARY: WATER
}
 

void PB_Parameters::initParams(double epssolvent, ofstream &os){


  this->epssolute = 1.0;
  this->tiny_real=0.0000000001;
  this->FPEPSI=gmath::physConst.get_four_pi_eps_i();
  this->PI=gmath::physConst.get_pi();
  this->tinfoil_boundary=0.0;
  this->vacuum_boundary=1.0;
  this->adjusted_boundary=epssolvent; 
  this->ddg_conv = -0.001;
  this->hradius = 0.05;
  this->nmax_surfpoints=15000;
  this->default_chshape = 0; //hat
  this->alpha1=0.15;
  this->alpha2=0.05;
  this->nalias1=2;
  this->nalias2=8;
  this->threadnum=1;
  this->FFTlambda=1.5;
  this->kappa=0.0; // kappa for Debye Hueckel boundary charge density in NPBC of FD
  this->debugvar=0; // for debugging: 0 or 1: no or yes
  this->convergence_fd=0.000001;;
  this->convergence_fft=0.3;

  // ewald edir default realcut,tolerance,kmax
  this->tole=0.000001;
  this->rcut=0.9;
  this->kmax=64;

  // quadrupole moment trace of spc water model
  this->quadr=0.0082;


  // xi ewald constant
  this->xiew=-2.837297;

  os << "# PB_PARAMS: epssolute " << epssolute << endl;
  os << "# PB_PARAMS: tiny_real " << tiny_real << endl;
  os << "# PB_PARAMS: tinfoil_boundary " << tinfoil_boundary << endl;
  os << "# PB_PARAMS: vacuum_boundary " << vacuum_boundary << endl;
  os << "# PB_PARAMS: adjusted_boundary " << adjusted_boundary << endl;
  os << "# PB_PARAMS: ddg_conv " << ddg_conv << endl;
  os << "# PB_PARAMS: hradius " << hradius << endl;
  os << "# PB_PARAMS: alpha1 " << alpha1 << endl;
  os << "# PB_PARAMS: alpha2 " << alpha2 << endl;
  os << "# PB_PARAMS: nalias1 " << nalias1 << endl;
  os << "# PB_PARAMS: nalias2 " << nalias2 << endl;
  os << "# PB_PARAMS: FFTlambda " << FFTlambda << endl;
  os << "# PB_PARAMS: kappa " << kappa << endl;
  os << "# PB_PARAMS: convergence_fd " << convergence_fd << endl;
  os << "# PB_PARAMS: convergence_fft " << convergence_fft << endl;
  os << "# PB_PARAMS: quadr " << quadr << endl;
  os << "# PB_PARAMS: tole " << tole << endl;
  os << "# PB_PARAMS: rcut " << rcut << endl;
  os << "# PB_PARAMS: kmax " << kmax << endl;
  os << "# PB_PARAMS: xiew " << xiew << endl;
}



double PB_Parameters::getPI(){return this->PI;}
double PB_Parameters::getFPEPSI(){return this->FPEPSI;};
double PB_Parameters::getEpssolute(){return this->epssolute;};
double PB_Parameters::getTiny_real(){return this->tiny_real;};
double PB_Parameters::get_rf_boundary_eps(){return this->adjusted_boundary;}; // adjusted boundary conditions
double PB_Parameters::get_sc_boundary_eps(){return this->vacuum_boundary;}; // vacuum boundary conditions
double PB_Parameters::get_ls_boundary_eps(){return this->tinfoil_boundary;}; // tinfoil boundary conditions
double PB_Parameters::get_hradius(){return this->hradius;};
int PB_Parameters::get_nmax_surfpoints(){return this->nmax_surfpoints;};
int PB_Parameters::get_default_chshape(){return this->default_chshape;};

double PB_Parameters::get_alpha1(){return this->alpha1;};
double PB_Parameters::get_alpha2(){return this->alpha2;};
int PB_Parameters::get_nalias1(){return this->nalias1;};
int PB_Parameters::get_nalias2(){return this->nalias2;};

double PB_Parameters::get_FFTlambda(){return this->FFTlambda;};

int PB_Parameters::get_threadnum(){return this->threadnum;};

double PB_Parameters::getKappa(){return this->kappa;};

int PB_Parameters::get_debugvar(){return this->debugvar;};

double PB_Parameters::get_convergence_fd(){return this->convergence_fd;};
double PB_Parameters::get_convergence_fft(){return this->convergence_fft;};


double PB_Parameters::get_rcut(){return this->rcut;};
double PB_Parameters::get_tole(){return this->tole;};
int PB_Parameters::get_kmax(){return this->kmax;};

double PB_Parameters::get_quadr(){return this->quadr;};
double PB_Parameters::get_xiew(){return this->xiew;};
