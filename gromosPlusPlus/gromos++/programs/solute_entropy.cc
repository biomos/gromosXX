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
 * @file solute_entropy.cc
 * Calculates Schlitter and Quasiharmonic analysis
 */
/**
 * @page programs Program Documentation
 *
 * @anchor solute_entropy
 * @section solute_entropy Schlitter and Quasiharmonic analysis
 * @author @ref rba @ref ns
 * @date 15-09-08
 *
 * Program solute_entropy takes a coordinate trajectory and calculates the
 * configurational entropy in @f$kJ\cdot K^{-1}\cdot mol^{-1}@f$ using the Schlitter and the quasiharmonic analysis 
 * methods for a given set of atoms.
 * The entropy can be averaged over a window of a given size.
 * If requested, a rotational fit prior to the entropy calculation is carried out.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td>[\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;]</td></tr>
 * <tr><td> \@atomsentropy</td><td>&lt;@ref AtomSpecifier "atoms" to consider for entropy&gt; </td></tr>
 * <tr><td> \@temp</td><td>&lt;temperature (K)&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * <tr><td>[\@ref</td><td>&lt;reference structure to fit against&gt;]</td></tr>
 * <tr><td>[\@ref_pbc</td><td>&lt;boundary type for reference for fit&gt;]</td></tr>
 * <tr><td>[\@atomsfit</td><td>&lt;@ref AtomSpecifier "atoms" to consider for fit&gt;]</td></tr>
 * <tr><td>[\@method</td><td>&lt;methods to use: schlitter quasiharm (default both)&gt;]</td></tr>
 * <tr><td>[\@n</td><td>&lt;entropy is calculated every nth step (default every step)&gt;]</td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  solute_entropy
    @topo       ex.top
    @pbc        r
    @time       0 0.1
    @atomsentropy  1:CB
    @atomsfit   1:CB
    @ref        ex.coo
    @ref_pbc    r
    @traj       ex.tr
    @temp       300
    @method     schlitter quasiharm
    @n          10
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <string>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector_double.h>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Physics.h"
#include "../src/utils/Rmsd.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/fit/PositionUtils.h"
#include "../src/utils/groTime.h"
#include "../src/gromos/Exception.h"

using namespace utils;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace gmath;
using namespace fit;
using namespace std;

// Function declaration ---------------------------------------------
// Function to compute the entropy ----------------------------------
double entropy(gsl_matrix * cov, gsl_vector * av, gsl_vector * mass, 
               unsigned int confs, int ndof, double temperature, unsigned int n_step);
// computing the entropy with the quasiharmonic analysis
double freq_etc(gsl_matrix * cov, gsl_vector * av, gsl_vector * mass, 
                unsigned int confs, int ndof, double temperature, unsigned int n_step);



int main(int argc, char **argv) {

  Argument_List knowns; 
  knowns << "topo" << "pbc" << "ref"  <<"ref_pbc" << "atomsfit" << "temp" 
         << "atomsentropy" << "n" << "time" << "method" << "traj"
          << "list";

  string usage = "# ";
  usage += string(argv[0]) + "\n";
  usage += "\n\t@topo        <topology>\n";
  usage += "\t@pbc           <boundary type>\n";
  usage += "\t[@list         <atom_list for gathering>]\n";
  usage += "\t@ref           <structure to fit against>\n";
  usage += "\t@ref_pbc       <boundary type for reference for fit>\n"; 
  usage += "\t@atomsfit      <atoms to consider for fit>\n";
  usage += "\t@atomsentropy  <atoms to consider for entropy>\n";
  usage += "\t@temp          <temperature (K)>\n";
  usage += "\t@time          <time and dt>\n";
  usage += "\t@method        <methods to use schlitter=yes, quasiharm=yes (default both)>\n";
  usage += "\t@n             <entropy calculated every nth step (default every step)>\n";
  usage += "\t@traj          <trajectory files>\n";


  try {
    
    // define necessary constants
    //Boltzmann constant in J/K
    const double KB = gmath::physConst.get_boltzmann() * 1000 / gmath::physConst.get_avogadro();
    const double E = gmath::physConst.get_euler(); /// Euler number
    const double HBAR = gmath::physConst.get_hbar() / (1e9 * gmath::physConst.get_avogadro()); /// Plank constant over 2 Pi, SI units
    const double MU = gmath::physConst.get_atomic_mass_unit(); /// Atomic mass unit
    const double NA = gmath::physConst.get_avogadro(); /// Avogadros number

    cout.precision(10);
    cout << "# " << KB << " " << E << " " << HBAR << " " << MU << " " << NA << endl;

    Arguments args(argc, argv, knowns, usage);


    //  read topology -------------------------------
    args.check("topo", 1);
    InTopology it(args["topo"]);
    System sys(it.system());

    System refSys(sys);

    // DW : if gathering based on an atom list, now read in the list
    Arguments::const_iterator pbciter = args.lower_bound("pbc");
    ++pbciter;

    string gath = pbciter->second;
    cout << "# gather option : " << gath << endl;

    if(gath=="1" || gath == "4"){
        if(args.count("list") == 0){
            cout << "################# Gathering : WARNING #################" << endl
                 << "  You have requested to gather the system based on an atom list," << endl
                 << "while you didn't define such a list, therefore the gathering"<< endl
                 << "will be done according to the 1st atom of the previous molecule." << endl;
        } else {
            AtomSpecifier gathlist(sys);

            if(args.count("list") > 0){
                Arguments::const_iterator iter = args.lower_bound("list");
                Arguments::const_iterator to = args.upper_bound("list");

                //int testid=0;
                for(;iter!=to;iter++){
                    string spec=iter->second.c_str();
                    gathlist.addSpecifierStrict(spec);
                }
                for(unsigned int j=0;j<gathlist.size()/2;++j){
                    int i=2*j;
                    sys.primlist[gathlist.mol(i)][0]=gathlist.atom(i);
                    sys.primlist[gathlist.mol(i)][1]=gathlist.mol(i+1);
                    sys.primlist[gathlist.mol(i)][2]=gathlist.atom(i+1);

                    refSys.primlist[gathlist.mol(i)][0]=gathlist.atom(i);
                    refSys.primlist[gathlist.mol(i)][1]=gathlist.mol(i+1);
                    refSys.primlist[gathlist.mol(i)][2]=gathlist.atom(i+1);

                    cout << "# updated prim : mol " << gathlist.mol(i) << " atom " << gathlist.atom(i)
                         << "# refe : mol " << sys.primlist[gathlist.mol(i)][1] << " atom " << sys.primlist[gathlist.mol(i)][2] << endl;

                }
            }
        }
    }
    // end here

    // parse the type of boundary conditions and create pbc
    Boundary *pbc = BoundaryParser::boundary(sys, args, "pbc");
    // parse the gathering method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args, "pbc");

    //   get simulation temperature ------------------
    double temp = args.getValue<double>("temp", true);

    //   get simulation time-step dt [ps]
    //   the trajectory is read every nread steps
    //   intervals for analysis: nana steps (0 -> only at the end)
    Time time(args);

    //   skip diagonalization for a number of steps?  ------------------
    unsigned int n_step = args.getValue<unsigned int>("n", false, 1);

    //  which method: schlitter, quasiharmonic, both -------------
    bool schlitter = true, quasi = true;
    {
      Arguments::const_iterator iter = args.lower_bound("method"),
      to = args.upper_bound("method");
      if (iter != to) {
        schlitter = false;
        quasi = false;
      }
      // search arguments for methods
      for(; iter != to; ++iter) {
        if (iter->second == "schlitter")
          schlitter = true;
        if (iter->second == "quasiharm")
          quasi = true;
        if (iter->second == "both")
          schlitter = quasi = true;
      }
    }
    if (!schlitter && !quasi) {
      throw gromos::Exception("solute_entropy", "@method You should do either schlitter or quasi.");
    }

    // only solute molecules -------------------------
    // use the reference class as container for the selected atoms
    Reference atoms_entropy(&sys);
    AtomSpecifier atoms(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("atomsentropy");
      Arguments::const_iterator to = args.upper_bound("atomsentropy");
      for (; iter != to; iter++) 
        atoms.addSpecifier(iter->second);
    }

    // add atoms to the reference
    atoms_entropy.addAtomSpecifier(atoms);
    if (atoms.empty()) 
      throw gromos::Exception("solute_entropy", "@atomsentropy no atoms specified.");

    int ndof = 3 * atoms.size(); // nr of degrees of freedom for entropy


    // fitting -------------------------------------------
    //System refsys(sys);
    Reference ref(&refSys);
    bool fit = false;
    {
      InG96 fitRef;
      Arguments::const_iterator iter = args.lower_bound("ref");
      if (iter != args.upper_bound("ref")) {
        fit = true;
        //InG96 fitRef;
        fitRef.open(iter->second.c_str());
        fitRef.select("SOLUTE");
      }else{
        //InG96 fitRef;
        fit = true;
        if(args.count("traj")>0)
            fitRef.open(args.lower_bound("traj")->second);
      }
        fitRef >> refSys;
        fitRef.close();
        // correct for periodic boundary conditions by calling
        // the appropriate gathering method	      
        // parse the type of boundary conditions and create pbc
        if (args.count("ref_pbc") > 0) {
          Boundary *refpbc = BoundaryParser::boundary(refSys, args, "ref_pbc");
          // parse the gathering method
          Boundary::MemPtr refgathmethod = args::GatherParser::parse(sys,refSys,args, "ref_pbc");
          (*refpbc.*refgathmethod)();
          delete refpbc;
        }
      //}
    }

    if (fit) {
      AtomSpecifier ref_atoms(refSys);
      Arguments::const_iterator iter = args.lower_bound("atomsfit");
      Arguments::const_iterator to = args.upper_bound("atomsfit");
      for (; iter != to; iter++)
        ref_atoms.addSpecifier(iter->second);

      //if (ref_atoms.empty())
      //  throw gromos::Exception("solute_entropy", "@atomsfit no atoms specified.");
      if (ref_atoms.empty()){
          cout << "# WARNING : @atomsfit not specified. Will use @atomsentropy for fit."
                  << endl;
          ref_atoms = atoms;
      }

      // add atoms to the reference
      ref.addAtomSpecifier(ref_atoms);
    }

    // define the necessary vectors and matrices ----------------
    // averages 
    gsl_vector *pos_av = gsl_vector_alloc(ndof);
    gsl_vector_set_zero(pos_av);
    // position vector to store coordinates intermediately
    gsl_vector *position = gsl_vector_alloc(ndof);
    // masses (seems to be 3 times too big 
    // -> leaves space for mass scaling of specific DOFs)
    gsl_vector *mass = gsl_vector_alloc(ndof);
    // matrices
    gsl_matrix * covariance = gsl_matrix_alloc(ndof, ndof);
    gsl_matrix_set_zero(covariance);


    // fill the mass vector
    int count = 0;
    for (int i = 0; i < sys.numMolecules(); ++i) {
      for (int j = 0; j < sys.mol(i).numAtoms(); ++j) {
        if (atoms_entropy.weight(i, j) > 0) {
          for (int k = 0; k < 3; ++k) {
            gsl_vector_set(mass, count, sys.mol(i).topology().atom(j).mass());
            count++;
          }
        }
      }
    }  
    
    cerr << "WARNING: this program is experimental." << endl;
    
    // define input coordinate
    InG96 ic;
    RotationalFit rf(&ref);
    unsigned int numFrames = 0;
    
    // loop over all trajectories
    for (Arguments::const_iterator iter = args.lower_bound("traj"),
            to = args.upper_bound("traj"); iter != to; ++iter) {

      ic.open((iter->second).c_str());

      while (!ic.eof()) {
        ic.select("SOLUTE");
        ic >> sys >> time;

        // correct for periodic boundary conditions by calling
        // the appropriate gathering method	 
        (*pbc.*gathmethod)();
        if (fit)
          rf.fit(&sys);

        // print title
        if (numFrames == 0) {
          cout.setf(ios::right, ios::adjustfield);
          cout << "# " << setw(15) << "time" << setw(15) << "schlitter"
               << setw(15) << "quasi" << endl;
        }
        
        // get the positions
        count = 0;
        for (int i = 0; i < sys.numMolecules(); ++i) {
          for (int j = 0; j < sys.mol(i).numAtoms(); ++j) {
            if (atoms_entropy.weight(i, j) > 0) {
              for (int k = 0; k < 3; ++k) {
                gsl_vector_set(position, count, sys.mol(i).pos(j)[k]);
                count++;
              } // k
            } // if selected
          } // atoms
        } // molecules


        // fill the average vector/matrix
        double temporary;
        for (int i = 0; i < ndof; i++) {
          temporary = gsl_vector_get(position, i);
          gsl_vector_set(pos_av, i, gsl_vector_get(pos_av, i) + temporary);
          for (int j = i; j < ndof; j++)
            gsl_matrix_set(covariance, i, j, (gsl_matrix_get(covariance, i, j) + temporary * gsl_vector_get(position, j)));
        }
        
        double entr_schlitter = 0.0, entr_quasi = 0.0;
        
        if (schlitter)
          entr_schlitter = entropy(covariance, pos_av, mass, numFrames, 
                                   ndof, temp, n_step);
        
        if (quasi)
          entr_quasi = freq_etc(covariance, pos_av, mass, numFrames, 
                                ndof, temp, n_step);
        if (entr_schlitter > 0.0 || entr_quasi > 0.0 || numFrames == 0) {
          cout.precision(2);
          cout.setf(ios::right, ios::adjustfield);
          cout << time;

          cout.precision(8);
          cout.setf(ios::fixed, ios::floatfield);
          cout << setw(15) << entr_schlitter << setw(15) << entr_quasi << endl;
        }
        
        numFrames++;
      } // frames
      ic.close();
    } // files

    gsl_matrix_free(covariance);
    gsl_vector_free(mass);
    gsl_vector_free(pos_av);
    gsl_vector_free(position);

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

// Function to compute the entropy ----------------------------------
double entropy(gsl_matrix * cov, gsl_vector * av, gsl_vector * mass, 
               unsigned int confs, int ndof, double temperature, unsigned int n_step) {

  // define necessary constants
  //Boltzmann constant in J/K
  const double KB = gmath::physConst.get_boltzmann() * 1000 / gmath::physConst.get_avogadro();
  const double E = gmath::physConst.get_euler(); /// Euler number
  const double HBAR = gmath::physConst.get_hbar() / (1e9 * gmath::physConst.get_avogadro()); /// Plank constant over 2 Pi, SI units
  const double MU = gmath::physConst.get_atomic_mass_unit(); /// Atomic mass unit
  const double NA = gmath::physConst.get_avogadro(); /// Avogadros number

  const double dbl_confs = double(confs+1);
  // calculate factor: mu*k*T*e^2/hbar^2 (+conversion from nm to m)
  double mukte2h2 = temperature * 1e-18 * MU * KB * E * E / (HBAR * HBAR);
  double JtoKJ = 1e-3;
  gsl_matrix * lu = gsl_matrix_alloc(ndof, ndof);

  double temp;
  // off-diagonal elements
  for (int i = 0; i < ndof; i++) {
    for (int j = i + 1; j < ndof; j++) {

      temp = gsl_matrix_get(cov, i, j) / dbl_confs - gsl_vector_get(av, i) * gsl_vector_get(av, j) / (dbl_confs * dbl_confs);
      temp = temp * mukte2h2;
      gsl_matrix_set(lu, i, j, temp * gsl_vector_get(mass, j));
      gsl_matrix_set(lu, j, i, temp * gsl_vector_get(mass, i));
    }
  }
  // diagonal elements
  for (int i = 0; i < ndof; i++) {
    temp = gsl_vector_get(av, i) / dbl_confs;
    temp = -temp * temp + gsl_matrix_get(cov, i, i) / dbl_confs;
    temp = 1 + temp * mukte2h2 * gsl_vector_get(mass, i);
    gsl_matrix_set(lu, i, i, temp);
  }

  double diag = 0.0;
  for (int i = 0; i < ndof; i++)
    diag += log(gsl_matrix_get(lu, i, i));
  // double entropy_d = 0.5 * KB * NA * diag;

  double entropy = 0.0;
  if (confs % n_step == 0) {
    gsl_permutation * p = gsl_permutation_alloc(ndof);
    int s;
    gsl_linalg_LU_decomp(lu, p, &s);
    double lndet = gsl_linalg_LU_lndet(lu);
    entropy = 0.5 * KB * NA * lndet;
    gsl_permutation_free(p);
  }
  
  // convert from J.K-1.mol-1 to kJ.K-1.mol-1
  entropy = entropy * JtoKJ;

  gsl_matrix_free(lu);

  return entropy;
}

double freq_etc(gsl_matrix * cov, gsl_vector * av, gsl_vector * mass, 
                unsigned int confs, int ndof, double temperature, unsigned int n_step) {

  // define necessary constants
  //Boltzmann constant in J/K
  const double KB = gmath::physConst.get_boltzmann() * 1000 / gmath::physConst.get_avogadro();
  //const double E = gmath::physConst.get_euler(); /// Euler number
  const double HBAR = gmath::physConst.get_hbar() / (1e9 * gmath::physConst.get_avogadro()); /// Plank constant over 2 Pi, SI units
  const double MU = gmath::physConst.get_atomic_mass_unit(); /// Atomic mass unit
  const double NA = gmath::physConst.get_avogadro(); /// Avogadros number

  const double dbl_confs = double(confs+1);
  
  // calculate factor: mu*k*T*e^2/hbar^2 (+conversion from nm to m)
  //double mukte2h2 = temperature * 1e-18 * MU * KB * E * E / (HBAR * HBAR);
  double mufac = 1e-18 * MU;
  double kT = temperature * KB;
  double JtoKJ = 1e-3;

  gsl_matrix * lu = gsl_matrix_alloc(ndof, ndof);

  double temp;
  // off-diagonal elements
  for (int i = 0; i < ndof; i++) {
    for (int j = i + 1; j < ndof; j++) {
      temp = gsl_matrix_get(cov, i, j) / dbl_confs - gsl_vector_get(av, i) * gsl_vector_get(av, j) / (dbl_confs * dbl_confs);
      temp = temp * mufac;
      // for a faster program -> use a sqrtmass_vector
      temp = temp * sqrt(gsl_vector_get(mass, i)) *
              sqrt(gsl_vector_get(mass, j));
      gsl_matrix_set(lu, i, j, temp);
      gsl_matrix_set(lu, j, i, temp);
    }
  }
  // diagonal elements
  for (int i = 0; i < ndof; i++) {
    temp = gsl_vector_get(av, i) / dbl_confs;
    temp = -temp * temp + gsl_matrix_get(cov, i, i) / dbl_confs;
    temp = temp * mufac * gsl_vector_get(mass, i);
    gsl_matrix_set(lu, i, i, temp);
  }

  double entropy = 0.0;
  if (confs % n_step == 0) {
    gsl_vector *eval = gsl_vector_alloc(ndof);
    gsl_matrix *evec = gsl_matrix_alloc(ndof, ndof);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(ndof);
    gsl_eigen_symmv(lu, eval, evec, w);
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

    double ho_kT;
    int nr_evals = 0;
    for (int i = 0; i < ndof; i++) {
      temp = sqrt(kT * gsl_vector_get(eval, i)) / HBAR;
      if (temp > 1.0e-10) {
        nr_evals++;
        ho_kT = 1.0 / temp;
        temp = KB * NA * (ho_kT / (exp(ho_kT) - 1.0) - log(1.0 - exp(-ho_kT)));
        entropy += temp;
      }
    }
    gsl_eigen_symmv_free(w);
    gsl_matrix_free(evec);
    gsl_vector_free(eval);
  }
  
  // convert from J.K-1.mol-1 to kJ.K-1.mol-1
  entropy = entropy * JtoKJ;
  
  gsl_matrix_free(lu);

  return entropy;
}

