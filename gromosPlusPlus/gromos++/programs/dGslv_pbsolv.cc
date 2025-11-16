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

/**
 * @file dGslv_pbsolv.cc
 * calculates the correction to the
 * electrostatic component of the
 * solvation free energy for either
 * going from LS/PBC to CB/NPBC or for
 * going from RF/PBC to CB/NPBC
 */

/**
 * @page programs Program Documentation
 *
 * @anchor dGslv_pbsolv
 * @section dGslv_pbsolv calculates a correction to the electrostatic component of the solvation free energy
 * @author @ref mr
 * @date 25-10-10
 *
 * Progam dGslv_pbsolv will compute two electrostatic components of the solvation free energy,
 * namely one for CB/NPBC and one for the user-specified electrostatics scheme (LS or RF).
 *
 * In the following, the abbreviation dGslv will be used to denote an electrostatic component of
 * the solvation free energy.
 *
 * dGslv will be computed for a user-specified group of atoms.
 * Note that all atoms of the solute topology block will be nonpolarizable,
 * i.e. they will be assigned a relative dielectric permittivity of one.
 *
 *
 * The CB/NPBC dGslv will be computed from a FD algorithm.
 *
 * 
 * The LS/PBC dGslv will be computed from both a FD and a FFT algorithm (for comparison).
 * But the user should use the FD value to cancel possible grid discretization errors.
 * That is, the resulting correction is: CB/NPBC[FD] - LS/PBC[FD]
 *
 *
 * The RF/PBC dGslv will be computed from a FFT algorithm.
 * But the user should use the corresponding LS calculation to compute a correction
 * to cancel possible grid discretization errors.
 * That is, the resulting correction is: CB/NPBC[FD] - LS/PBC[FD] + LS/PBC[FFT] - RF/PBC[FFT]
 *
 *
 * In the LS-scheme, tinfoil boundary conditions are used and a hat charge shaping function
 * will be used.
 *
 * In the RF-scheme, a user-specified relative dielectric permittivity is used.
 * Note that a relative dielectric permittivity of one implies no application of a reaction-fied correction.
 *
 *
 *
 * The solute will be centered in the computational box, with its center of geometry.
 *
 * The algorithms employed are from these papers:
 * FD: Comput. Phys. Commun. 62, 187-197 (1991)
 * FFT: J. Chem. Phys. 116, 7434-7451 (2002), J. Chem. Phys. 119, 12205-12223 (2003)
 *
 *
 *  *
 * Example:
 * @verbatim
 dGslv_pbsolv 
 @topo topo.top
 @atoms 1:a
 @atomsTOcharge 1:a
 @coord coor.dat
 @schemeELEC RF
 @epsSOLV  66.6
 @epsRF 66.6
 @epsNPBC 78.4
 @rcut 1.4
 @gridspacing 0.02
 @gridpointsXYZ 158 158 158
 @maxiter 600
 @cubesFFT 4
 @probeIAC 5
 @probeRAD 0.14
 @HRAD 0.05
 @pbc r gbond
 @radscal 1.0
 @rminORsigma 0
 @endverbatim
 *
 *
 *
 */


#ifdef HAVE_LIBFFTW3
print asmsadf $%;
#include <cassert>
#include <cctype>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <string>

#include <fftw3.h>
#include <bits/stl_bvector.h>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/AtomicRadii.h"
#include "../src/fit/PositionUtils.h"
#include "../src/pb/FDPoissonBoltzmann.h"
#include "../src/pb/FDPoissonBoltzmann_ICCG_NPBC.h"
#include "../src/pb/FDPoissonBoltzmann_ICCG_PBC.h"
#include "../src/pb/PB_Parameters.h"
#include "../src/pb/FFTPoisson.h"
#include "../src/pb/FFTGridType.h"
#include "../src/pb/FFTBoundaryCondition.h"
#include "../src/gio/InPDB.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gromos/Exception.h"

#include "../config.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace gmath;
using namespace bound;
using namespace args;
using namespace utils;
using namespace pb;

vector <double> fd_ls_npbc_slv(utils::AtomSpecifier atoms, utils::AtomSpecifier atomsTOcharge, int ngrid_x, int ngrid_y, int ngrid_z, double gridspacing, double epsNPBC,  int maxiter, double convergence_fd, double &result_npbc_slv, double &gridcenterX, double &gridcenterY, double &gridcenterZ, double &gridstartX, double &gridstartY, double &gridstartZ, ofstream &os) {

  // Create vectors for storing the calculated potentials
  vector <double> potentials_npbc_slv (0);

  FDPoissonBoltzmann_ICCG_NPBC iccg_npbc(ngrid_x,ngrid_y,ngrid_z);
  os << "# ************************************************** " <<   endl;  
  os << "# *** FD LS NPBC EPSSOLV *** " <<   endl;
  os << "# ************************************************** " <<   endl;
  
  FDPoissonBoltzmann pbsolv_NPBC_epssolvent(atoms,atomsTOcharge,ngrid_x,ngrid_y,ngrid_z, \
					    gridspacing, false,epsNPBC, os);
  pbsolv_NPBC_epssolvent.setupGrid(true, os);
  pbsolv_NPBC_epssolvent.solveforpotential_npbc(maxiter, convergence_fd,iccg_npbc, os);
  result_npbc_slv = pbsolv_NPBC_epssolvent.dGelec(os, &potentials_npbc_slv);
  return potentials_npbc_slv;
    
}

vector <double> fd_ls_npbc_vac (utils::AtomSpecifier atoms, utils::AtomSpecifier atomsTOcharge, int ngrid_x, int ngrid_y, int ngrid_z, double gridspacing, double epssolvent, int maxiter, double convergence_fd, double &result_npbc_vac, double &gridcenterX, double &gridcenterY, double &gridcenterZ, double &gridstartX, double &gridstartY, double &gridstartZ, ofstream &os){
  vector <double> potentials_npbc_vac (0);
  FDPoissonBoltzmann_ICCG_NPBC iccg_npbc(ngrid_x,ngrid_y,ngrid_z);
  
  os << "# ************************************************** " <<   endl;  
  os << "# *** FD LS NPBC VAC *** " <<   endl;
  os << "# ************************************************** " <<   endl;
  
  FDPoissonBoltzmann pbsolv_NPBC_vac(atoms,atomsTOcharge,ngrid_x,ngrid_y,ngrid_z, \
				     gridspacing, false,1.0, os);
  pbsolv_NPBC_vac.setupGrid(true, os, gridstartX, gridstartY, gridstartZ, gridcenterX, gridcenterY, gridcenterZ);
  pbsolv_NPBC_vac.solveforpotential_npbc(maxiter, convergence_fd,iccg_npbc, os);
  result_npbc_vac = pbsolv_NPBC_vac.dGelec(os, &potentials_npbc_vac);
 
  return potentials_npbc_vac;
}

vector <double> fd_ls_pbc_slv(utils::AtomSpecifier atoms, utils::AtomSpecifier atomsTOcharge, int ngrid_x, int ngrid_y, int ngrid_z, double gridspacing, double epssolvent, int maxiter, double convergence_fd, double &result_pbc_slv, double &gridcenterX, double &gridcenterY, double &gridcenterZ, double &gridstartX, double &gridstartY, double &gridstartZ, ofstream &os) {
  vector <double> potentials_pbc_slv (0);
  FDPoissonBoltzmann_ICCG_PBC iccg_pbc(ngrid_x,ngrid_y,ngrid_z);
  
  os << "# ************************************************** " <<   endl;
  os << "# *** FD LS PBC EPSSOLV *** " <<   endl;
  os << "# ************************************************** " <<   endl;
  
  
  FDPoissonBoltzmann pbsolv_PBC_epssolvent(atoms,atomsTOcharge,ngrid_x,ngrid_y,ngrid_z,\
					   gridspacing, true,epssolvent, os);
  pbsolv_PBC_epssolvent.setupGrid(true, os, gridstartX, gridstartY, gridstartZ, gridcenterX, gridcenterY, gridcenterZ);
  pbsolv_PBC_epssolvent.solveforpotential_pbc(maxiter, convergence_fd,iccg_pbc, os);
  result_pbc_slv =  pbsolv_PBC_epssolvent.dGelec(os, &potentials_pbc_slv);
  os << "# ************************************************** " <<   endl;
  
  return potentials_pbc_slv;
}

vector <double> fd_ls_pbc_vac(utils::AtomSpecifier atoms, utils::AtomSpecifier atomsTOcharge,	\
			      int ngrid_x, int ngrid_y, int ngrid_z, double gridspacing,\
			      int maxiter, double convergence_fd, double &result_pbc_vac, double &gridcenterX, double &gridcenterY, double &gridcenterZ, double &gridstartX, double &gridstartY, double &gridstartZ, 
			      			\
			      ofstream &os) {
  vector <double> potentials_pbc_vac (0);
  FDPoissonBoltzmann_ICCG_PBC iccg_pbc(ngrid_x,ngrid_y,ngrid_z);
  
  os << "# ************************************************** " <<   endl;
  os << "# *** FD LS PBC VAC *** " <<   endl;
  os << "# ************************************************** " <<   endl;
  
  FDPoissonBoltzmann pbsolv_PBC_vac(atoms,atomsTOcharge,ngrid_x,ngrid_y,ngrid_z,gridspacing, true,1.0, os); //vacuum calc
  pbsolv_PBC_vac.setupGrid(true, os, gridstartX, gridstartY, gridstartZ, gridcenterX, gridcenterY, gridcenterZ);
  pbsolv_PBC_vac.solveforpotential_pbc(maxiter, convergence_fd,iccg_pbc, os);
  result_pbc_vac = pbsolv_PBC_vac.dGelec(os, &potentials_pbc_vac);

  return potentials_pbc_vac;
}

vector <double> fft_ls_pbc(utils::AtomSpecifier atoms, utils::AtomSpecifier atomsTOcharge, int ngrid_x, int ngrid_y, int ngrid_z, double gridspacing, double epssolvent, PB_Parameters ppp, double rcut, double epsRF, int maxiter, double convergence_fft, int fftcub, ofstream &os) {
  vector <double> potentials_fft_ls_pbc (0);

  FFTGridType gt(ngrid_x, ngrid_y, ngrid_z,				\
		 ngrid_x*gridspacing, ngrid_y*gridspacing, ngrid_z*gridspacing, fftcub, os);
  
  os << "# ************************************************** " <<   endl;
  os << "# *** FFT LS (PBC) *** " <<   endl;
  os << "# ************************************************** " <<   endl;
  
  
  // make the boundary conditions for LS
  FFTBoundaryCondition bc_LS(0, "LS",
			     ppp.get_alpha1(), ppp.get_alpha2(), ppp.get_nalias1(), ppp.get_nalias2(), rcut, epsRF, os);
  // (rcut and epsRF are basically dummies)
  // and print them
  bc_LS.dumpparameters(os);
  
  // setup the main object
  FFTPoisson fftp_LS(atoms, atomsTOcharge, gt, bc_LS, gridspacing, maxiter, convergence_fft, ppp.get_FFTlambda(),
		     epssolvent, false, true, os);
  
  // and now we iterate
  os << "# call solve_poisson ..." << endl;
  fftp_LS.solve_poisson(os, &potentials_fft_ls_pbc);
  
  return potentials_fft_ls_pbc;
}

    //*****************************************************************
vector <double> fft_rf_pbc(utils::AtomSpecifier atoms, utils::AtomSpecifier atomsTOcharge, int ngrid_x, int ngrid_y, int ngrid_z, double gridspacing, double epssolvent, PB_Parameters ppp, double rcut, double epsRF, int maxiter, double convergence_fft, int fftcub, ofstream &os) {
  vector <double> potentials_fft_rf_pbc (0);

  // DO AN ADDITIONAL RF FFT
  FFTGridType gt(ngrid_x, ngrid_y, ngrid_z,				\
		 ngrid_x*gridspacing, ngrid_y*gridspacing, ngrid_z*gridspacing, fftcub, os);
  
  os << "# ************************************************** " <<   endl;
  os << "# *** FFT RF (PBC) *** " <<   endl;
  os << "# ************************************************** " <<   endl;
  
  if (fabs(epsRF- 1.0)> ppp.tiny_real ){
    os << "# WITH REACTION FIELD CORRECTION : RF ... "  << endl;
    // make the boundary conditions for RF
    FFTBoundaryCondition bc_RF(1, "RF",
			       ppp.get_alpha1(), ppp.get_alpha2(), ppp.get_nalias1(), ppp.get_nalias2(), rcut, epsRF, os);
    // setup the main objects for RF
    FFTPoisson fftp_RF(atoms, atomsTOcharge,  gt, bc_RF, gridspacing, maxiter, convergence_fft, ppp.get_FFTlambda(), epssolvent, false, false, os);
    // print params and iterate : RF
    bc_RF.dumpparameters(os);
    os << "# call solve_poisson ..." << endl;
    fftp_RF.solve_poisson(os, &potentials_fft_rf_pbc);  
  }
  else{
    os << "# WITHOUT REACTION FIELD CORRECTION : SC ... "  << endl;
    // make the boundary conditions for SC
    FFTBoundaryCondition bc_SC(2, "SC",
			       ppp.get_alpha1(), ppp.get_alpha2(), ppp.get_nalias1(), ppp.get_nalias2(), rcut, epsRF, os);
    // setup the main objects for RF
    FFTPoisson fftp_SC(atoms,  atomsTOcharge,  gt, bc_SC, gridspacing, maxiter, convergence_fft, ppp.get_FFTlambda(),
		       epssolvent, false, false, os);
    // print params and iterate : SC
    bc_SC.dumpparameters(os);
    os << "# call solve_poisson ..." << endl;
    fftp_SC.solve_poisson(os);
  }
  
  return potentials_fft_rf_pbc;
  
} 
   
void writeout(string schemeELEC, utils::AtomSpecifier atomsTOcharge, vector <double> potentials_npbc_slv, vector <double> potentials_npbc_vac, vector <double> potentials_pbc_slv, vector <double> potentials_pbc_vac, vector <double> potentials_fft_ls_pbc, vector <double> potentials_fft_rf_pbc) {
  // WRITEOUT
  cout << setw(6) << "# atom" << setw(7) << "nres" << setw(7) << "name" << setw(10) << "charge" << setw(16) << "NPBC_SLV" << setw(16) << "NPBC_VAC" << setw(16) << "PBC_SLV" << setw(16) << "PBC_VAC";
  if (schemeELEC == "RF"){
    cout << setw(16) << "FFT_LS_PBC" << setw(16) << "FFT_RF_PBC";
    }
  cout << endl;
  
  for (unsigned i = 0; i < atomsTOcharge.size(); ++i) {
    cout << setw(6) << i+1
	 << setw(7) << atomsTOcharge.resnum(i)+1
	 << setw(7) << atomsTOcharge.name(i)
	 << setw(10) << fixed << std::setprecision(4) << atomsTOcharge.charge(i)
	 << setw(16) << fixed << std::setprecision(7) << potentials_npbc_slv.at(i)
	 << setw(16) << fixed << std::setprecision(7) << potentials_npbc_vac.at(i)
	 << setw(16) << fixed << std::setprecision(7) << potentials_pbc_slv.at(i)
	 << setw(16) << fixed << std::setprecision(7) << potentials_pbc_vac.at(i);
    if (schemeELEC == "RF"){
      cout  << setw(16) << fixed << std::setprecision(7) << potentials_fft_ls_pbc.at(i)
	    << setw(16) << fixed << std::setprecision(7) << potentials_fft_rf_pbc.at(i);
    }
    cout << endl;
  }
  cout << endl;
  // END WRITEOUT
  }


int main(int argc, char **argv){
  
  Argument_List knowns;
  knowns << "topo" 
         << "pbc"
         << "atoms" <<  "atomsTOcharge" << "coord" << "pqr" << "schemeELEC" << "epsSOLV"
         << "epsRF" << "rcut"
         << "gridspacing" << "coordinates" << "maxiter" << "nogridpoints" << "NPBCsize"
         << "cubesFFT" << "probeIAC" << "probeRAD" << "HRAD" <<  "epsNPBC" <<  "radscal" << "rminORsigma" << "increasegrid" << "verbose";

  string usage = "# " + string(argv[0]);
  usage += "\n\n# USAGE\n";
  usage += "\n";
  usage += "# -----------------------------------------------------------------------------------------\n";
  usage += "# if you have a gromos coordinate file:\n";
  usage += "\t@topo            <molecular topology file>\n";
  usage += "\t@pbc             <boundary type> [<gather method>]\n";
  usage += "\t@coord           <gromos96 coordinates>\n";
  usage += "\t@probeIAC        <integer atom code for radius calculation\n";
  usage += "\t                  (for water oxygen it would be 4 or 5, depending on the ff)>\n";
  usage += "\t@atoms           <atoms to include for the pb calculations\n";
  usage += "\t                  expected in gromos format>\n";
  usage += "\t@atomsTOcharge   <atoms to charge; expected in gromos format>\n";
  usage += "\t@rminORsigma     <how to calculate the radii - rmin (0) or sigma (1); default: 0>\n";
  usage += "\n";
  usage += "# -----------------------------------------------------------------------------------------\n";
  usage += "# if you have a pqr file:\n";
  usage += "# short description - the last two elements in each ATOM or HETATOM line\n";
  usage += "# are the atom charge and the atom radii (in Angstrom);\n";
  usage += "# hydrogen atoms can have zero radius (see @radH below)\n";
  usage += "\t@pqr             <pqr file>\n";
  usage += "\t@coordinates     <box coordinates in X Y Z direction (in nm) that were used\n";
  usage += "\t                  in the simulation>\n";
  usage += "\t@atoms           <atoms to include for the pb calculations>\n";
  usage += "\t                  typically, all atoms should be included\n";
  usage += "\t                  (simply type 'a' to include all atoms)\n";
  usage += "\t                  if not, atoms are indexed from 1 and individual atoms\n";
  usage += "\t                  can be picked by\n";
  usage += "\t                  using a comma delimiter; ranges can be specified using a hyphen\n";
  usage += "\t                  e.g. '1,2,3,4,5,8,9' is the same as '1-5,8,9'>\n";
  usage += "\t@atomsTOcharge   <atoms to charge; e.g. the perturbed atoms in the\n";
  usage += "\t                  free energy calculations;\n";
  usage += "\t                  same rules as above to pick individual atoms>\n";
  usage += "\n";
  usage += "# -----------------------------------------------------------------------------------------\n";
  usage += "# general input:\n";
  usage += "\t@schemeELEC      <electrostatics scheme: LS or RF>\n";
  usage += "\t if RF: @rcut    <cutoff distance in nm>\n";
  usage += "\t if RF: @epsRF   <dielectric permittivity of the reaction field>\n";
  usage += "\t@epsSOLV         <relative dielectric permittivity of the employed solvent model>\n";

  usage += "\t@gridspacing     <grid spacing in nm; should not be much higher than 0.02 but\n";
  usage += "\t                  lower numbers are computationally much more expensive in\n";
  usage += "\t                  terms of time and memory usage>\n";
  usage += "\t[@nogridpoints    <optional; if number is given gridspacing is ignored and the\n";
  usage += "\t                  grid is set up based on the specified number of grids, where\n";
  usage += "\t                  the number of gridpoints is the same in each direction x,y,z;\n";
  usage += "\t                  default the number of gridpoints is calculated based on the box-size\n";
  usage += "\t                  and the gridspacing>]\n";
  usage += "\t[@NPBCsize        <optional; if argument is given change the NPBC grid size, based\n";
  usage += "\t                  on the box dimension given in nm; default the npbc box is\n";
  usage += "\t                  extended by 4 nm>]\n";
  usage += "\t[@epsNPBC        <relative dielectrict permittivity for calculation of macroscopic,\n";
  usage += "\t                  non-periodic boundary conditions; default 78.4 (for water)>]\n";
  usage += "\t[@maxiter        <maximum number of iteration steps; default 600>]\n";
  usage += "\t[@cubesFFT       <number of cubes in the fast Fourier transformation for\n";
  usage += "\t                  boundary smoothing; default 4>]\n";
  usage += "\t[@probeRAD       <probe radius in nm; default 0.14 (for water)>]\n";
  usage += "\t[@radH           <your desired hydrogen radius in nm; default: 0.05>]\n";

  usage += "\t[@radscal        <scale non-H radii with this factor (use this only in case you\n";
  usage += "\t                  want to play with radii); default 1.0>]\n";
  usage += "\t[@increasegrid   <takes three integer values for X Y Z; grid for PBC calculations gets increased by the number of given gridpoints;\n";
  usage += "\t                  may be usefull if atoms close to the border of the box extend the grid!>]\n";
  usage += "\t[@verbose        <path to log file to document status and errors>]\n";
  
  try{

    Arguments args(argc, argv, knowns, usage);
    
    // ------------------------------------
    // READ NON-SYSTEM DEPENDENT PARAMETERS
    // ------------------------------------
    // turn verbose mode on
    std::ofstream os;
    bool verbose=0;
    if(args.count("verbose")>=0) {
      string log_file = "";
      log_file = args["verbose"].c_str();
      if (log_file=="")  throw gromos::Exception("dGslv_pbsolv","verbose - no file name given");
      os.open(log_file.c_str());
    } else {
      os.open("/dev/null");
    }
    
    
    // read schemeELEC
    string schemeELEC = "";
    if(args.count("schemeELEC")>0){
      schemeELEC = args["schemeELEC"];
      transform(schemeELEC.begin(), schemeELEC.end(), schemeELEC.begin(), static_cast<int (*)(int)>(std::toupper));
    }
    if(schemeELEC!="LS" && schemeELEC !="RF")
      throw gromos::Exception("dGslv_pbsolv","schemeELEC format "+schemeELEC+" unknown. Exiting ...");  
    os << "# READ: schemeELEC " << schemeELEC << endl;

    // read cutoff distance
    double rcut=-1;
    if(schemeELEC=="RF") {
      if(args.count("rcut")>0) rcut=atof(args["rcut"].c_str());
      if (rcut<0)  throw gromos::Exception("dGslv_pbsolv","No or negative RF cutoff given. Exiting ...");
      os << "# READ: rcut " << rcut << endl;
    }

    // read probe radius
    double probe_rad=0.14;
    if(args.count("probeRAD")>0) probe_rad=atof(args["probeRAD"].c_str());
    if (probe_rad<0)  throw gromos::Exception("dGslv_pbsolv","The probe radius (probeRAD) must not be negative. Exiting ...");
    os << "# READ: probe_rad " << probe_rad << endl;

    // read hydrogen radius
    double hydrogen_rad=0.05;
    if(args.count("radH")>0) hydrogen_rad=atof(args["radH"].c_str());
    if (hydrogen_rad<0)  throw gromos::Exception("dGslv_pbsolv","The hydrogen radius (radH) must not be negative. Exiting ...");
    os << "# READ: hydrogen_rad " << hydrogen_rad << endl;

    // read radscal
    double radscal=1.0;
    if(args.count("radscal")>0) radscal=atof(args["radscal"].c_str());
    os << "# READ: radscal " << radscal << endl;

    // read increasegrid
    int increasegrid_x=0;
    int increasegrid_y=0;
    int increasegrid_z=0;
    if(args.count("increasegrid")>0) {
      Arguments::const_iterator iterincreasegrid=args.lower_bound("increasegrid");
      if(iterincreasegrid!=args.upper_bound("coordinates")){
	increasegrid_x=atoi(iterincreasegrid->second.c_str());
      ++iterincreasegrid;
      }
      if(iterincreasegrid!=args.upper_bound("coordinates")){
	increasegrid_y=atoi(iterincreasegrid->second.c_str());
	++iterincreasegrid;
      }
      if(iterincreasegrid!=args.upper_bound("coordinates")){
	increasegrid_z=atoi(iterincreasegrid->second.c_str());
      }
      os << "# READ: increasegrid " << increasegrid_x << " " << increasegrid_y << " " << increasegrid_z << endl;
    }
    // read epsilon
    double epssolvent=0.0;
    if(args.count("epsSOLV")>0) epssolvent=atof(args["epsSOLV"].c_str());
    if (epssolvent<1)  throw gromos::Exception("dGslv_pbsolv","The solvent permittivity (epsSOLV) not given or given value is smaller than 1. Exiting ...");
    os << "# READ: epssolvent " << epssolvent << endl;

    // read RF epsilon
    double epsRF=-1;
    if(schemeELEC=="RF") {
      if(args.count("epsRF")>0) epsRF=atof(args["epsRF"].c_str());
      if (epsRF<1 && epsRF != 0.0)  throw gromos::Exception("dGslv_pbsolv","Reaction field permittivity (epsRF) not given or given value not allowed. Exiting ...");
      os << "# READ: epsRF " << epsRF << endl;
    }

    // read NPBC epsilon
    double epsNPBC=78.4;
    if(args.count("epsNPBC")>0) epsNPBC=atof(args["epsNPBC"].c_str());
    if (epsNPBC<1)  throw gromos::Exception("dGslv_pbsolv","The permittivity for the non-periodic boundary conditions (epsNPBC) must not be smaller than 1. Exiting ...");
    os << "# READ: epsNPBC " << epsNPBC << endl;

    // read maxiter
    int maxiter=600;
    if(args.count("maxiter")>0) maxiter=atoi(args["maxiter"].c_str());
    if (maxiter<=0)  throw gromos::Exception("dGslv_pbsolv","The maximum number of iterations (maxiter) must be positive. Exiting ...");
    os << "# READ: maxiter " << maxiter << endl;

    // read radius definition
    int rminorsigma=0;
    if(args.count("coord")>0) {
      if(args.count("rminORsigma")>0) rminorsigma=atoi(args["rminORsigma"].c_str());
      if (rminorsigma != 0 && rminorsigma != 1)  throw gromos::Exception("dGslv_pbsolv","The rminORsigma flag should be 0 (rmin) or 1 (sigma). Exiting ...");
      os << "# READ: rminORsigma " << rminorsigma << endl;
    }

    // read cubesFFT
    int fftcub=4;
    if(args.count("cubesFFT")>0) fftcub=atoi(args["cubesFFT"].c_str());
    if (fftcub<=0)  throw gromos::Exception("dGslv_pbsolv","The FFT cubelet number must be positive. Exiting ...");
    os << "# READ: fftcub " << fftcub << endl;

    // ------------------------------------------------
    // FINISHED READING NON-SYSTEM DEPENDENT PARAMETERS
    // ------------------------------------------------

    // -----------------
    // CASE CNF IS GIVEN
    // -----------------
    
    if(args.count("coord")>0) {
    // read topology
    args.check("topo",1);
    gio::InTopology cnf_in_top(args["topo"]);
    System cnf_sys(cnf_in_top.system());  // here the system is created based on the topology
    System cnf_refSys(cnf_in_top.system());

    // set atoms for which we want to compute the solvation free energy
    utils::AtomSpecifier cnf_atoms(cnf_sys);
    utils::AtomSpecifier cnf_atomsTOcharge(cnf_sys);

    Arguments::const_iterator iter1=args.lower_bound("atoms");
    Arguments::const_iterator to=args.upper_bound("atoms");
    for(;iter1!=to;iter1++){
      cnf_atoms.addSpecifier(iter1->second.c_str());
    } 
     
    if (cnf_atoms.size()==0)  throw gromos::Exception("dGslv_pbsolv","No atoms specified. Exiting ...");

    iter1=args.lower_bound("atomsTOcharge");
    to=args.upper_bound("atomsTOcharge");
    for(;iter1!=to;iter1++){
      cnf_atomsTOcharge.addSpecifier(iter1->second.c_str());
    }
  

    if (cnf_atomsTOcharge.size()==0)  throw gromos::Exception("dGslv_pbsolv","No atoms to charge specified. Exiting ...");

    os << "# READ: atomsTOcharge " << endl;
    for (unsigned int i=0;i< cnf_atomsTOcharge.size();i++){
      os << "# READ: mol " <<  cnf_atomsTOcharge.mol(i) << " " << cnf_atomsTOcharge.name(i) << endl;
    }

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(cnf_sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(cnf_sys,cnf_refSys,args);


    // read  coordinates
      
    InG96 ic(args["coord"]);
    if(args.count("coord")>0){
      ic.select("ALL"); //it reads in ALL molecules (including solvent)
      ic >> cnf_sys; // that's the input stream
      (*pbc.*gathmethod)();
      ic.close();
    }
    if (!cnf_sys.hasPos)throw gromos::Exception("dGslv_pbsolv","No position block in coordinate file. Exiting ...");


    // calculate number of gridpoints either with argument @nogridpoints or @gridspacing
    double a=cnf_sys.box().K()[0];
    double b=cnf_sys.box().L()[1];
    double c=cnf_sys.box().M()[2];
    if (!cnf_sys.hasBox) throw gromos::Exception("dGslv_pbsolv","No box block in coordinate file. Exiting ...");
    if (a <= 0) throw gromos::Exception("dGslv_pbsolv","At least one coordinate is zero. Exiting ...");
    if (b <= 0) throw gromos::Exception("dGslv_pbsolv","At least one coordinate is zero. Exiting ...");
    if (c <= 0) throw gromos::Exception("dGslv_pbsolv","At least one coordinate is zero. Exiting ...");
  
    int ngrid_x=0;
    int ngrid_y=0;
    int ngrid_z=0;
    
    // read nogridpoints
    int ngrid=0;
    if(args.count("nogridpoints")>0) {
      ngrid=atof(args["nogridpoints"].c_str());
      if (ngrid<0 && ngrid != 0)  throw gromos::Exception("dGslv_pbsolv","The number of grids must not be 0 or negative. Exiting ...");
      os << "# READ: nogridpoints " << ngrid << endl;
      ngrid_x=ngrid;
      ngrid_y=ngrid;
      ngrid_z=ngrid;
    }
    
    // read gridspacing
    double gridspacing=0.0;
    if(args.count("gridspacing")>0) {
      gridspacing=atof(args["gridspacing"].c_str());
      if (gridspacing<0)  throw gromos::Exception("dGslv_pbsolv","The grid spacing must not be negative. Exiting ...");
      os << "# READ: gridspacing " << gridspacing << endl;
      // calculate grid size for periodic boxes
      ngrid_x=ceil(a/gridspacing);
      ngrid_y=ceil(b/gridspacing);
      ngrid_z=ceil(c/gridspacing);
    }

    // adapt gridspacing such that it fits periodic box or calculate gridspacing if @nogridpoints argument is given
    //gridspacing= a/ngrid_x;


    // calculate grid size for non-periodic boxes
    // we use the boxdimensions+4nm or the the specified boxsize in nm given by the argument NPBCsize -
    // dependent on if NPBCsize is given
    int ngrid_x_npbc=0;
    int ngrid_y_npbc=0;
    int ngrid_z_npbc=0;
    double NPBCsize=0.0;

    if(args.count("NPBCsize")>0) {
      NPBCsize=atof(args["NPBCsize"].c_str());
      if (NPBCsize<0 && NPBCsize != 0)  throw gromos::Exception("dGslv_pbsolv","The box size for NPBC must not be 0 or negative. Exiting ...");
      os << "# READ: NPBCsize " << NPBCsize << endl;
      ngrid_x_npbc=ceil((NPBCsize)/gridspacing);
      ngrid_y_npbc=ceil((NPBCsize)/gridspacing);
      ngrid_z_npbc=ceil((NPBCsize)/gridspacing);
    }
    else {
      ngrid_x_npbc=ceil((a+4)/gridspacing); // or ngrid_x_npbc=ceil(ngrid_x+(4/gridspacing)-1);
      ngrid_y_npbc=ceil((b+4)/gridspacing); // or ngrid_y_npbc=ceil(ngrid_y+(4/gridspacing)-1);
      ngrid_z_npbc=ceil((c+4)/gridspacing); // or ngrid_z_npbc=ceil(ngrid_z+(4/gridspacing)-1);
    }



    // scale the box (upon user specification)
    if (increasegrid_x > 0) {
      ngrid_x += increasegrid_x;
      os << "Increased the grid on the x axis by " << increasegrid_x << " according to @increasegrid flag!" << endl;
    }
    if (increasegrid_y > 0) {
      ngrid_y += increasegrid_y;
      os << "Increased the grid on the y axis by " << increasegrid_y << " according to @increasegrid flag!" << endl;
    }
    if (increasegrid_z > 0) {
      ngrid_z += increasegrid_z;
      os << "Increased the grid on the z axis by " << increasegrid_z << " according to @increasegrid flag!" << endl;
    }

    os << "# Boxsize: " << a << endl;
    os << "# Adapted gridspacing: " << gridspacing << endl;
    os << "# Calculated number of gridpoints PBC (x,y,z): " << ngrid_x << " " << ngrid_y << " " << ngrid_z << endl;
    os << "# Calculated number of gridpoints NPBC (x,y,z): " << ngrid_x_npbc << " " << ngrid_y_npbc << " " << ngrid_z_npbc << endl;


    // set PB parameters
    PB_Parameters ppp(epssolvent, os);
    // get convergence
    double convergence_fd=ppp.get_convergence_fd();
    double convergence_fft=ppp.get_convergence_fft();

    // read probe IAC
    int probe_iac=0;
    if(args.count("probeIAC")>0) probe_iac=atoi(args["probeIAC"].c_str());
    if (probe_iac<=0)  throw gromos::Exception("dGslv_pbsolv","The probe integer atom code (probeIAC) must not be negative or 0. Exiting ...");
    os << "# READ: probe_iac " << probe_iac << endl;

    
    // get atomic radii (these are rmin)
    os << "# RADII: based on probe_iac " << probe_iac << endl;
    if (rminorsigma == 1 ){
      os << "# RADII: based on sigma" << endl;
    }
    else
      {
	os << "# RADII: based on rmin" << endl;
      }
    utils::compute_atomic_radii_vdw(probe_iac-1,probe_rad,cnf_sys, cnf_in_top.forceField());


    // loop over atoms and set H radii to hradius_param
    // also, if sigma is used, adapt appropriately
    for (unsigned int i=0;i<cnf_atoms.size();i++){
  
      if (fabs(cnf_atoms.radius(i))< ppp.tiny_real){
	cnf_sys.mol(cnf_atoms.mol(i)).topology().atom(cnf_atoms.atom(i)).setradius( hydrogen_rad );
	//os << "# !!!! H atom number (0 radius) " << i << " : assign rad(i) = " << cnf_atoms.radius(i) << endl;
      }
      else{
	// scale radius
	if ( rminorsigma == 1){
	  // use sigma instead (i.e. add probe_rad again, then divide by 2^(1/6) then subtract probe_rad)
	  double tmprad=(cnf_atoms.radius(i) + probe_rad)/( exp((1.0/6.0)  * log(2.0) ) )-probe_rad;
	  // now scale
	  tmprad=tmprad*radscal;
	  cnf_sys.mol(cnf_atoms.mol(i)).topology().atom(cnf_atoms.atom(i)).setradius( tmprad);
	}
	else{
	  // rmin ok, just scale
	  double tmprad=cnf_atoms.radius(i)*radscal;
	  cnf_sys.mol(cnf_atoms.mol(i)).topology().atom(cnf_atoms.atom(i)).setradius( tmprad);
	}
      }
      //os << "# radius of atom " << i << " : rad(i) = " << cnf_atoms.radius(i)  << endl;
    }
    for (unsigned int i=0;i<cnf_atomsTOcharge.size();i++){
      //os << "# radius of atom_to_charge " << i << " : rad(i) = " << cnf_atomsTOcharge.radius(i)  << endl;
    }

    //start the machine
      double result_npbc_slv = 0;
      double result_npbc_vac = 0;
      double result_pbc_slv = 0;
      double result_pbc_vac = 0;
      double gridcenterX = cnf_sys.box().K()[0]/2;
      double gridcenterY = cnf_sys.box().L()[1]/2;
      double gridcenterZ = cnf_sys.box().M()[2]/2;
      double gridstartX = 0;
      double gridstartY = 0;
      double gridstartZ = 0;
      
      vector <double> potentials_npbc_slv;
      vector <double> potentials_npbc_vac;
      vector <double> potentials_pbc_slv;
      vector <double> potentials_pbc_vac;
      vector <double> potentials_fft_ls_pbc;
      vector <double> potentials_fft_rf_pbc;
      
      potentials_npbc_vac = fd_ls_npbc_vac(cnf_atoms, cnf_atomsTOcharge, ngrid_x_npbc, ngrid_y_npbc, ngrid_z_npbc, gridspacing, epssolvent, maxiter, convergence_fd, result_npbc_vac, gridcenterX, gridcenterY, gridcenterZ, gridstartX, gridstartY, gridstartZ, os);
      potentials_npbc_slv = fd_ls_npbc_slv(cnf_atoms, cnf_atomsTOcharge, ngrid_x_npbc, ngrid_y_npbc, ngrid_z_npbc, gridspacing, epsNPBC, maxiter, convergence_fd, result_npbc_slv, gridcenterX, gridcenterY, gridcenterZ, gridstartX, gridstartY, gridstartZ, os);
      double result_ls_npbc = result_npbc_slv - result_npbc_vac;
      os << "# DGRESULT NPBC " << result_ls_npbc << endl;

      gridstartX = 0;
      gridstartY = 0;
      gridstartZ = 0;
      
      potentials_pbc_vac = fd_ls_pbc_vac(cnf_atoms, cnf_atomsTOcharge, ngrid_x, ngrid_y, ngrid_z, gridspacing, maxiter, convergence_fd, result_pbc_vac, gridcenterX, gridcenterY, gridcenterZ, gridstartX, gridstartY, gridstartZ, os);
      potentials_pbc_slv = fd_ls_pbc_slv(cnf_atoms, cnf_atomsTOcharge, ngrid_x, ngrid_y, ngrid_z, gridspacing, epssolvent, maxiter, convergence_fd, result_pbc_slv, gridcenterX, gridcenterY, gridcenterZ, gridstartX, gridstartY, gridstartZ, os);
      double result_ls_pbc = result_pbc_slv - result_pbc_vac;
      os << "# DGRESULT PBC " << result_ls_pbc << endl;
      
      if (schemeELEC == "RF" ) {
      potentials_fft_ls_pbc = fft_ls_pbc(cnf_atoms, cnf_atomsTOcharge, ngrid_x, ngrid_y, ngrid_z, gridspacing, epssolvent, ppp, rcut, epsRF, maxiter, convergence_fft, fftcub, os);

      potentials_fft_rf_pbc = fft_rf_pbc(cnf_atoms, cnf_atomsTOcharge, ngrid_x, ngrid_y, ngrid_z, gridspacing, epssolvent, ppp, rcut, epsRF, maxiter, convergence_fft, fftcub, os);
      }
      
      writeout(schemeELEC, cnf_atomsTOcharge, potentials_npbc_slv, potentials_npbc_vac, potentials_pbc_slv, potentials_pbc_vac, potentials_fft_ls_pbc, potentials_fft_rf_pbc);


    
    os.close();




    }
    // -----------------
    // CASE PQR IS GIVEN
    // -----------------
    // when pqr is given, a molecule topology of argon atoms is created
    // these argon atoms have the charges and the radii from the pqr file
    // then, the system is created from this molecule topology
    // also, a "dummy-solvent" molecule is added to the system
    // because the system needs a solvent molecule (but it's not needed for any calculations)
    else if (args.count("pqr")>0){

      // --------------------
      // read in the pqr file
      InPDB pqr_in(args["pqr"]);
      pqr_in.select("ALL");
      pqr_in.readPQR();

      // -------------------------------------
      // create an atom topology of an Ar atom
      AtomTopology pqr_atom_top; 
      pqr_atom_top.setIac(31); // we generate a fake Ar atom


      // ---------------------------------------------------------------------------------------
      // create a molecule topology; assign atom topologies with charges and radii from pqr file
      MoleculeTopology pqr_mol_top;
      for (unsigned int i = 0; i < pqr_in.numAtoms(); i++) {
	pqr_atom_top.setName(pqr_in.getAtomName(i));
	pqr_atom_top.setCharge(pqr_in.PQR_getCharges(i));
	pqr_atom_top.setradius(pqr_in.PQR_getRadii(i)/10.0); // div by 10 from A to nm
	pqr_mol_top.addAtom(pqr_atom_top);
      }

      // --------------------------------
      // create a "fake" solvent topology
      SolventTopology pqr_solv_top;
      pqr_solv_top.addAtom(pqr_atom_top);

      // ----------------------------------------
      // create a Molecule and assign coordinates
      Molecule pqr_mol(pqr_mol_top);
      pqr_mol.initPos();
      for (unsigned int i = 0; i < pqr_in.numAtoms(); i++) {
	pqr_mol.pos(i) = pqr_in.getAtomPos(i)/10.0; // div by 10 for conv from A in nm
      }

      // ----------------------------------------------
      // create a "fake" solvent (needs no coordinates)
      Solvent pqr_solv(pqr_solv_top);
      
      // -------------------------------------------------------------
      // Finally, create a system and add the molecule and the solvent
      System pqr_sys;
      pqr_sys.addMolecule(pqr_mol);
      pqr_sys.addSolvent(pqr_solv);

      // -------------------------------------------------------------------------------
      // set atoms used and atoms for which we want to compute the solvation free energy
      // if pqr is read in, user is allowed to type atoms in gromos format (mol:atoms)
      // but also just atoms are accepted (e.g. 1-40);
      // if molecule is not specified, '1:' is added to the string
      // (as there is only one molecule when read from a pqr)
      utils::AtomSpecifier pqr_atoms(pqr_sys);
      utils::AtomSpecifier pqr_atomsTOcharge(pqr_sys);

      string atoms_to_add;
      // do it for the atoms
      Arguments::const_iterator iter1=args.lower_bound("atoms");
      Arguments::const_iterator to=args.upper_bound("atoms");
      for(;iter1!=to;iter1++){
	atoms_to_add = iter1->second.c_str();
	if (atoms_to_add.find(":") == -1) { // find returns -1 if char not found
	  atoms_to_add = "1:" + atoms_to_add;
	}
	pqr_atoms.addSpecifier(atoms_to_add);
      }
      if (pqr_atoms.size()==0)  {
	throw gromos::Exception("dGslv_pbsolv","No atoms specified. Exiting ...");
      }
      // do it for the atomsTOcharge
      iter1=args.lower_bound("atomsTOcharge");
      to=args.upper_bound("atomsTOcharge");
      for(;iter1!=to;iter1++){
	atoms_to_add = iter1->second.c_str();
	if (atoms_to_add.find(":") == -1) { // find returns -1 if char not found
	  atoms_to_add = "1:" + atoms_to_add;
	}
	pqr_atomsTOcharge.addSpecifier(atoms_to_add);
      }
      
      if (pqr_atomsTOcharge.size()==0) {
	throw gromos::Exception("dGslv_pbsolv","No atoms to charge specified. Exiting ...");
      }
      
      os << "# READ: atomsTOcharge " << endl;
      for (unsigned int i=0;i< pqr_atomsTOcharge.size();i++){
      	os << "# READ: mol " <<  pqr_atomsTOcharge.mol(i) << " " << pqr_atomsTOcharge.name(i) << endl;
      }

      
      // ----------------
      // read coordinates and calculate the number of gridpoints needed
      double a=0.0;
      double b=0.0;
      double c=0.0;
      int ngrid_x=0;
      int ngrid_y=0;
      int ngrid_z=0;
      
      Arguments::const_iterator iter2=args.lower_bound("coordinates");
      if(iter2!=args.upper_bound("coordinates")){
	a=atof(iter2->second.c_str());
      ++iter2;
      }
      if(iter2!=args.upper_bound("coordinates")){
	b=atof(iter2->second.c_str());
	++iter2;
      }
      if(iter2!=args.upper_bound("coordinates")){
	c=atof(iter2->second.c_str());
      }
      if (a<=0 || b<=0 || c<=0) throw gromos::Exception("dGslv_pbsolv","Coordinates must be positive and given in X Y Z (space-separated, no comma inbetween). Exiting ...");

      

    if (a <= 0) throw gromos::Exception("dGslv_pbsolv","At least one coordinate is zero. Exiting ...");
    if (b <= 0) throw gromos::Exception("dGslv_pbsolv","At least one coordinate is zero. Exiting ...");
    if (c <= 0) throw gromos::Exception("dGslv_pbsolv","At least one coordinate is zero. Exiting ...");

    int ngrid=0;
    
    // read nogridpoints
    if(args.count("nogridpoints")>0) {ngrid=atof(args["nogridpoints"].c_str());
      if (ngrid<0 && ngrid != 0)  throw gromos::Exception("dGslv_pbsolv","The number of grids must not be 0 or negative. Exiting ...");
      os << "# READ: nogridpoints " << ngrid << endl;
      ngrid_x=ngrid;
      ngrid_y=ngrid;
      ngrid_z=ngrid;
    }
    
    // read gridspacing
    double gridspacing=0.0;
    if(args.count("gridspacing")>0) {
      gridspacing=atof(args["gridspacing"].c_str());
      if (gridspacing<0)  throw gromos::Exception("dGslv_pbsolv","The grid spacing must not be negative. Exiting ...");
      os << "# READ: gridspacing " << gridspacing << endl;

      // calculate grid size for periodic boxes
      ngrid_x=ceil(a/gridspacing);
      ngrid_y=ceil(b/gridspacing);
      ngrid_z=ceil(c/gridspacing);
    }

    // adapt gridspacing such that it fits periodic box or calculate gridspacing if @nogridpoints argument is given
    gridspacing= a/ngrid;

    // calculate grid size for non-periodic boxes
    // we use the boxdimensions+4nm or the the specified boxsize in nm given by the argument NPBCsize -
    // dependent on if NPBCsize is given
    int ngrid_x_npbc=0;
    int ngrid_y_npbc=0;
    int ngrid_z_npbc=0;
    double NPBCsize=0.0;

    if(args.count("NPBCsize")>0) {
      NPBCsize=atof(args["NPBCsize"].c_str());
      if (NPBCsize<0 && NPBCsize != 0)  throw gromos::Exception("dGslv_pbsolv","The box size for NPBC must not be 0 or negative. Exiting ...");
      os << "# READ: NPBCsize " << NPBCsize << endl;
      ngrid_x_npbc=ceil((NPBCsize)/gridspacing);
      ngrid_y_npbc=ceil((NPBCsize)/gridspacing);
      ngrid_z_npbc=ceil((NPBCsize)/gridspacing);

    }
    else {
      ngrid_x_npbc=ceil((a+4)/gridspacing); // or ngrid_x_npbc=ceil(ngrid_x+(4/gridspacing)-1);
      ngrid_y_npbc=ceil((b+4)/gridspacing); // or ngrid_y_npbc=ceil(ngrid_y+(4/gridspacing)-1);
      ngrid_z_npbc=ceil((c+4)/gridspacing); // or ngrid_z_npbc=ceil(ngrid_z+(4/gridspacing)-1);
    }
      
      // scale the box (upon user specification)

      if (increasegrid_x > 0) {
	ngrid_x += increasegrid_x;
	os << "Increased the grid on the x axis by " << increasegrid_x << " according to @increasegrid flag!" << endl;
      }
      if (increasegrid_y > 0) {
	ngrid_y += increasegrid_y;
	os << "Increased the grid on the y axis by " << increasegrid_y << " according to @increasegrid flag!" << endl;
      }
      if (increasegrid_z > 0) {
	ngrid_z += increasegrid_z;
	os << "Increased the grid on the z axis by " << increasegrid_z << " according to @increasegrid flag!" << endl;
      }
      
    os << "# Calculated number of gridpoints PBC (x,y,z): " << ngrid_x << " " << ngrid_y << " " << ngrid_z << endl;
    os << "# Calculated number of gridpoints NPBC (x,y,z): " << ngrid_x_npbc << " " << ngrid_y_npbc << " " << ngrid_z_npbc << endl;
    
      
      // -----------------
      // set PB parameters
      PB_Parameters ppp(epssolvent, os);
      // get convergence
      double convergence_fd=ppp.get_convergence_fd();
      double convergence_fft=ppp.get_convergence_fft();


      // ------------------------------------------------
      // loop over atoms and set H radii to hradius_param
      // also, if sigma is used, adapt appropriately
      for (unsigned int i=0;i<pqr_atoms.size();i++){
	
	if (fabs(pqr_atoms.radius(i))< ppp.tiny_real){
	  pqr_sys.mol(pqr_atoms.mol(i)).topology().atom(pqr_atoms.atom(i)).setradius( hydrogen_rad );
	  os << "# !!!! H atom number (0 radius) " << i << " : assign rad(i) = " \
	     << pqr_atoms.radius(i) << endl;
	}
	else{
	  // scale radius
	  if ( rminorsigma == 1){
	    // use sigma instead
	    // (i.e. add probe_rad again, then divide by 2^(1/6) then subtract probe_rad)
	    double tmprad=(pqr_atoms.radius(i) + probe_rad)/( exp((1.0/6.0)  * log(2.0) ) )-probe_rad;
	    // now scale
	    tmprad=tmprad*radscal;
	    pqr_sys.mol(pqr_atoms.mol(i)).topology().atom(pqr_atoms.atom(i)).setradius( tmprad);
	  }
	  else{
	    // rmin ok, just scale
	    double tmprad=pqr_atoms.radius(i)*radscal;
	    pqr_sys.mol(pqr_atoms.mol(i)).topology().atom(pqr_atoms.atom(i)).setradius( tmprad);
	  }
	}
	os << "# radius of atom " << i << " : rad(i) = " << pqr_atoms.radius(i)  << endl;
      }
      for (unsigned int i=0;i<pqr_atomsTOcharge.size();i++){
	os << "# radius of atom_to_charge " << i << " : rad(i) = " \
	   << pqr_atomsTOcharge.radius(i)  << endl;
      }


      // -----------------
      // start the machine

      double result_npbc_slv = 0;
      double result_npbc_vac = 0;
      double result_pbc_slv = 0;
      double result_pbc_vac = 0;
      double gridcenterX = a/2;
      double gridcenterY = b/2;
      double gridcenterZ = c/2;
      double gridstartX = 0;
      double gridstartY = 0;
      double gridstartZ = 0;
      
      vector <double> potentials_npbc_slv;
      vector <double> potentials_npbc_vac;
      vector <double> potentials_pbc_slv;
      vector <double> potentials_pbc_vac;
      vector <double> potentials_fft_ls_pbc;
      vector <double> potentials_fft_rf_pbc;
      
      potentials_npbc_slv = fd_ls_npbc_slv(pqr_atoms, pqr_atomsTOcharge, ngrid_x_npbc, ngrid_y_npbc, ngrid_z_npbc, gridspacing, epsNPBC, maxiter, convergence_fd, result_npbc_slv, gridcenterX, gridcenterY, gridcenterZ, gridstartX, gridstartY, gridstartZ, os);
      
      potentials_npbc_vac = fd_ls_npbc_vac(pqr_atoms, pqr_atomsTOcharge, ngrid_x_npbc, ngrid_y_npbc, ngrid_z_npbc, gridspacing, epssolvent, maxiter, convergence_fd, result_npbc_vac, gridcenterX, gridcenterY, gridcenterZ, gridstartX, gridstartY, gridstartZ, os);

      double result_ls_npbc = result_npbc_slv - result_npbc_vac;
      os << "# DGRESULT NPBC " << result_ls_npbc << endl;

      gridstartX = 0;
      gridstartY = 0;
      gridstartZ = 0;

      
      potentials_pbc_slv = fd_ls_pbc_slv(pqr_atoms, pqr_atomsTOcharge, ngrid_x, ngrid_y, ngrid_z, gridspacing, epssolvent, maxiter, convergence_fd, result_pbc_slv, gridcenterX, gridcenterY, gridcenterZ, gridstartX, gridstartY, gridstartZ, os);
      
      potentials_pbc_vac = fd_ls_pbc_vac(pqr_atoms, pqr_atomsTOcharge, ngrid_x, ngrid_y, ngrid_z, gridspacing, maxiter, convergence_fd, result_pbc_vac, gridcenterX, gridcenterY, gridcenterZ, gridstartX, gridstartY, gridstartZ, os);

      double result_ls_pbc = result_pbc_slv - result_pbc_vac;
      os << "# DGRESULT PBC " << result_ls_pbc << endl;
      
      if (schemeELEC == "RF" ) {
      potentials_fft_ls_pbc = fft_ls_pbc(pqr_atoms, pqr_atomsTOcharge, ngrid_x, ngrid_y, ngrid_z, gridspacing, epssolvent, ppp, rcut, epsRF, maxiter, convergence_fft, fftcub, os);

      potentials_fft_rf_pbc = fft_rf_pbc(pqr_atoms, pqr_atomsTOcharge, ngrid_x, ngrid_y, ngrid_z, gridspacing, epssolvent, ppp, rcut, epsRF, maxiter, convergence_fft, fftcub, os);
      }
      
      writeout(schemeELEC, pqr_atomsTOcharge, potentials_npbc_slv, potentials_npbc_vac, potentials_pbc_slv, potentials_pbc_vac, potentials_fft_ls_pbc, potentials_fft_rf_pbc);
      
      os.close();
       
    
      std::exit(1);
        }
    else {
      throw gromos::Exception("dGslv_pbsolv","No input file. Use @coord or @pqr");
      }
    
  }   // end of try for argument reading
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
   }
}

// no FFTW
#else

#include <iostream>

int main()
{
  std::cerr << "\nconfigure could not find the FFTW libraries" << std::endl
	    << "needed to run this program." << std::endl << std::endl
	    << "You need to add them to your CPPFLAGS, CXXFLAGS, LDFLAGS" << std::endl
            << "or run ./configure --with-fftw=<path>" << std::endl << std::endl
	    << "Reconfigure and recompile to use this program" << std::endl;
  return 1;
}

#endif
