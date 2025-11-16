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
 * @file hvap.cc
 * Calculates the heat of vaporistaion from a coordinate trajectory
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor hvap
 * @section hval calculates the heat of vaporistaion from a coordinate trajectory
 * @author @ref ae
 * @date 07.07.2011
 *
 * Program hvap ...
 * 
 * NOTE: the current program version does only work for solute atoms/molecules.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@method</td><td>&lt;method to goarse grain: atomic or molecular&gt; </td></tr>
 * <tr><td> \@dist</td><td>&lt;min max ngrid&gt; </td></tr>
 * <tr><td> \@beads</td><td>&lt;number of atoms per bead (atomic)&gt; or </td></tr>
 * <tr><td>        </td><td>&lt;sequence of bead size within one molecule (molecular)&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; &lt;gather method&gt; </td></tr>
 * <tr><td> \@trc</td><td>&lt;simulation trajectory or coordinate file&gt;</td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
   rdf
     @topo   ex.top
 * @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_vector.h>

#include "../src/args/Arguments.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gcore/Box.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/LJType.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gmath/Physics.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/System.h"
#include "../src/utils/groTime.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gromos/Exception.h"

using namespace args;
using namespace bound;
using namespace gcore;
using namespace gio;
using namespace gmath;
using namespace std;
using namespace utils;

double LJ(int iac1, int iac2, double r2, GromosForceField &gff);
double Coulomb(double c1, double c2, double r2, double eps, double kap, double cut);
void putIntoBox(Vec &v, Box &box);

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "rf" << "T" << "trc" << "time";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo     <molecular topology file>\n";
  usage += "\t[@pbc     <boundary type (read from GENBOX block if not specified)> [<gather method>]>]\n";
  usage += "\t@rf       <cut off radiuse> <epsilon> <kappa>\n";
  usage += "\t[@time    <time and dt, used to overwrite the time of read from the trajectory file>]\n";
  usage += "\t@T        <temperature, used for the RT volume expansion term>\n";
  usage += "\t@trc      <simulation trajectory or coordinate file>\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // handle the time
    Time time(args);
    
    // read the temperature
    if(args.count("T") != 1) {
      throw gromos::Exception(argv[0], "check @T argument to specify the temperature");
    }
    double T;
    {
      stringstream ss;
      ss << args["T"];
      ss >> T;
      if(ss.fail() || ss.bad()) {
        stringstream msg;
        msg << "could not convert " << args["T"] << " to be used as temperature";
        throw gromos::Exception(argv[0], msg.str());
      }
    }
    
    // read reaction fiel parameters
    double eps, kap, cut;
    if (args.count("rf") != 3) {
      //throw gromos::Exception(argv[0], "too few arguments for the reaction-field parameters (@rf)");
    } else {
      Arguments::const_iterator it = args.lower_bound("rf");
      stringstream ss;
      ss << it->second;
      ss >> cut;
      ++it;
      ss << it->second;
      ss >> eps;
      ++it;
      ss << it->second;
      ss >> kap;
    }

    // read topology
    args.check("topo", 1);
    InTopology it(args["topo"]);
    System sys(it.system());
    GromosForceField gff = it.forceField();

    if (args.count("trc") < 1) {
      throw gromos::Exception(argv[0], "no coordinate or trajectory file specified (@trc)");
    }

    vector<double> hvaps;

    // print header
    cout << "#" << setw(19) << "time" << setw(20) << "H_vap" << endl;
    
    Arguments::const_iterator trcfirs = args.lower_bound("trc");
    Arguments::const_iterator trclast = args.upper_bound("trc");
    for (args::Arguments::const_iterator trc = trcfirs;
            trc != trclast; ++trc) {

      // the input coordinates
      InG96 ic;

      // the boundary
      bound::Boundary *pbc;

      // read boundary type, either from @pbc or GENBOX block
      if (args.count("pbc") > 0) { // read from arguments
        pbc = args::BoundaryParser::boundary(sys, args);
      } else { // read from GENBOX block
        if (args.count("pbc") == 0) {
          cerr << "WARNING: @pbc given with no argument(s), reading boundary type "
                  "from GENBOX block of the trajectory/coordinate file\n";
        }
        ic.open(trc->second.c_str());
        ic.select("SOLUTE");
        ic >> sys;
        pbc = args::BoundaryParser::boundary(sys);
        ic.close();
      }

      // loop over the configurations of the trajectory file
      ic.open(trc->second.c_str());
      ic.select("ALL");
      int numMol = sys.numMolecules();
      while (!ic.eof()) {

        double hvap = 0.0;

        ic >> sys >> time;

        // some numnbers needed in case tere is solvent
        int totAtSolv = sys.sol(0).numAtoms();
        int numAtSolv = sys.sol(0).topology().numAtoms();
        int numSolvMol = totAtSolv / numAtSolv;

        // loop over the intermolecular atom pairs
        int numAt1, numAt2;
        int a1, a2;
        int iac1, iac2;
        Vec pos1, pos2;
        double charge1, charge2;
        int m2;
        double r2;
        int s;
#ifdef OMP
#pragma omp parallel for private(numAt1, numAt2, a1, a2, iac1, iac2, pos1, pos2, charge1, charge2, m2, r2, s) reduction(+:hvap) schedule(dynamic)
#endif
        for (int m1 = 0; m1 < numMol; ++m1) {
          numAt1 = sys.mol(m1).numAtoms();
          for (a1 = 0; a1 < numAt1; ++a1) {
            iac1 = sys.mol(m1).topology().atom(a1).iac();
            pos1 = sys.mol(m1).pos(a1);
            charge1 = sys.mol(m1).topology().atom(a1).charge();
            // the solute-solute unteractions
            for (m2 = m1 + 1; m2 < numMol; ++m2) {
              numAt2 = sys.mol(m2).numAtoms();
              for (a2 = 0; a2 < numAt2; ++a2) {
                iac2 = sys.mol(m2).topology().atom(a2).iac();
                pos2 = sys.mol(m2).pos(a2);
                r2 = (pos1 - pbc->nearestImage(pos1, pos2, sys.box())).abs2();
                charge2 = sys.mol(m2).topology().atom(a2).charge();
                hvap += LJ(iac1, iac2, r2, gff);
                hvap += Coulomb(charge1, charge2, r2, eps, kap, cut);
              }
            }
            // the solute-solvent interactions
            for (s = 0; s < totAtSolv; ++s) {
              iac2 = sys.sol(0).topology().atom(s % numAtSolv).iac();
              pos2 = sys.sol(0).pos(s);
              r2 = (pos1 - pbc->nearestImage(pos1, pos2, sys.box())).abs2();
              charge2 = sys.sol(0).topology().atom(s % numAtSolv).charge();
              hvap += LJ(iac1, iac2, r2, gff);
              hvap += Coulomb(charge1, charge2, r2, eps, kap, cut);
            }
          }
        }
        // and the solvent-solvent interactions
        int s1, s2;
#ifdef OMP
#pragma omp parallel for private(s1, s2, iac1, iac2, pos1, pos2, r2, charge1, charge2) reduction(+:hvap) schedule(dynamic)
#endif
        for (int s = 0; s < totAtSolv; s += numAtSolv) {
          for (s1 = s; s1 < s + numAtSolv; ++s1) {
            for (s2 = s + numAtSolv; s2 < totAtSolv; ++s2) {
              iac1 = sys.sol(0).topology().atom(s1 % numAtSolv).iac();
              iac2 = sys.sol(0).topology().atom(s2 % numAtSolv).iac();
              pos1 = sys.sol(0).pos(s1);
              pos2 = sys.sol(0).pos(s2);
              r2 = (pos1 - pbc->nearestImage(pos1, pos2, sys.box())).abs2();
              charge1 = sys.sol(0).topology().atom(s1 % numAtSolv).charge();
              charge2 = sys.sol(0).topology().atom(s2 % numAtSolv).charge();
              hvap += LJ(iac1, iac2, r2, gff);
              hvap += Coulomb(charge1, charge2, r2, eps, kap, cut);
            }
          }
        }

        // divide by the number of molecules (solute and solvent)
        hvap /= (numMol + numSolvMol) + physConst.get_boltzmann() * T;
        hvaps.push_back(hvap);
        
        // write the time series
        cout.precision(9);
        cout << fixed << setw(20) << time.time() << scientific << setw(20) << hvap << endl;

      } // end of loop over configuration of the current trajectory file
    }

    double hvap = 0.0;
    for (unsigned int i = 0; i < hvaps.size(); ++i) {
      hvap += hvaps[i];
    }
    hvap /= hvaps.size();

    cout << "# average: " << hvap << endl;

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

double LJ(int iac1, int iac2, double r2, GromosForceField &gff) {
  double r6 = r2 * r2 * r2;
  AtomPair ap(iac1, iac2);
  double c12 = gff.ljType(ap).c12();
  double c6 = gff.ljType(ap).c6();
  return (c6 - c12 / r6) / r6;
}

double Coulomb(double c1, double c2, double r2, double eps, double kap, double cut) {
  if (c1 == 0 || c2 == 0) {
    return 0.0;
  } else {
    double Crf = (2 - 2 * eps)*(1 + kap * cut) - eps * (kap * cut)*(kap * cut) /
            ((1 + 2 * eps)*(1 + kap * cut) + eps * (kap * cut)* (kap * cut));
    double r = sqrt(r2);
    double cut3 = cut * cut * cut;
    return -(c1 * c2 * physConst.get_four_pi_eps_i() * ((1 / r) - (0.5 * Crf * r2 / cut3) - (1 - 0.5 * Crf / cut)));
  }
}
