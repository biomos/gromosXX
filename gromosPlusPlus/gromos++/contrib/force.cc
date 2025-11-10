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
 * @file force.cc
 * Calculates the total electrostatic force on atoms of group A from atoms of group B. 
 * 
 * @page contrib Contrib program documentation
 *
 * @anchor force
 * @section force calculates the LJ and CRF force between two groups of atoms
 * @author @ref ae
 * @date January 21, 2013
 * 
 * Program force calculates the total force (Lennard-Jones and Coulomb plus reaction field)
 * of atoms within group B on atoms of group A. If group A consists of more than one
 * atom the printed force is averaged over all atoms within group A.
 * 
 * The output of the program list the LJ, CRF, and total force of each configuration of the
 * input trajectories as a vector (@f$x@f$-, @f$y@f$-, and @f$z@f$-component). In addition, the projection
 * of theses force vectors on a specified vector (\@projvec) are printed if requested.
 * 
 * Two algorithms to calculate the pair list of atoms from group A to atoms of group B may
 * be chosen: CHARGEGROUP and ATOMIC. Note that this selection only affects the atoms
 * of group B since the pair list is always calculated from one single atom of group
 * A to a single or a group of atoms within group B.
 * 
 * \@cutmin allows for a minimum cut-off distance > 0.0, i.e. to only consider interactions
 * within a distance @f$r@f$ with
 * @f[ 0 <= cutmin <= r < cut @f]
 * where @f$cutmin@f$ is set to 0.0 if not specified otherwise (\@cutmin).
 * different.
 * 
 * This program is parallelised using OpenMP.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> [\@pbc</td><td>&lt;periodic boundary conditions&gt;]</td></tr>
 * <tr><td> \@pairlist</td><td>&lt;type (CHARGEGROUP or ATOMIC)&gt; </td></tr>
 * <tr><td> \@cut</td><td>&lt;cut-off radius for the force calculations&gt; </td></tr>
 * <tr><td> [\@cutmin</td><td>&lt;only consider pairs i,j with @f$r_ij >= cutmin@f$&gt;]</td></tr>
 * <tr><td> \@epskap</td><td>&lt;epsilon and kappa to be used Coulomb/RF calculation&gt; </td></tr>
 * <tr><td> \@atomsA</td><td>&lt;atoms of group A (atom specifier)&gt; </td></tr>
 * <tr><td> \@atomsB</td><td>&lt;atoms of group B (atom specifier)&gt; </td></tr>
 * <tr><td> [\@projvec</td><td>&lt;a vector specifier to project the force vector to (vector specifier)&gt;]</td></tr>
 * <tr><td> \@verbose</td><td>&lt;(prints some information about the variables used)&gt; </td></tr>
 * <tr><td> \@trc</td><td>&lt;positional simulation trajectory files&gt; </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
 force
    @topo     ex.top
    @pbc      r
    @pairlist CHARGEGROUP
    @cut      1.4
    @cutmin   1.0
    @epskap   78.5 0.0
    @atomsA   1:1
    @atomsB   1:a not(1:1)
    @projvec  atom(3:4,14)
    @trc      ex.trc.gz
 @endverbatim
 * 
 * <hr>
 */

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/System.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/bound/Boundary.h"
#include "../src/args/GatherParser.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/utils/groTime.h"
#include "../src/utils/SimplePairlist.h"
#include "../src/gmath/Physics.h"
#include "../src/utils/VectorSpecifier.h"
#include "../src/utils/Value.h"
#include "../src/gcore/LJType.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace args;
using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace gmath;

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "pairlist" << "cut" << "cutmin" << "cutrf" << "epskap"
          << "atomsA" << "atomsB" << "trc" << "verbose" << "projvec";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t[@pbc     <periodic boundary conditions>]\n"
          "\t          (only useful in case you want to overrule the entries of GENBOX)\n";
  usage += "\t@pairlist <type (CHARGEGROUP or ATOMIC)>\n";
  usage += "\t@cut      <cut-off radius for the force calculations>\n";
  usage += "\t[@cutmin  <only consider pairs i,j with r_ij >= cutmin>]\n";
  usage += "\t[@cutrf   <to set the cutoff of the reaction field different from cut>]\n";
  usage += "\t@epskap   <epsilon and kappa to be used Coulomb/RF calculation>\n";
  usage += "\t@atomsA   <atoms of group A (atom specifier)>\n";
  usage += "\t@atomsB   <atoms of group B (atom specifier)>\n";
  usage += "\t[@projvec  <a vector specifier to project the force vector to (vector specifier)>]\n";
  usage += "\t           (e.g. atom(1:1,2) the vector pointing from atom 2 to 1,\n";
  usage += "\t                 cart(x,y,z) the vector with Cartesian coordinates x, y, and z)\n";
  usage += "\t@verbose  (prints some information about the variables used)\n";
  usage += "\t@trc      <positional simulation trajectory files>\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology
    args.check("topo", 1);
    InTopology it(args["topo"]);
    System sys(it.system());
    System refSys(it.system()); // a reference topology, in our case the same as the actual topology

    // The GROMOS force field as read from the topology
    GromosForceField gff(it.forceField());

    // parse boundary conditions (from GENBOX, if no @pbc is given, from @pbc else)
    Boundary *pbc = BoundaryParser::boundary(sys, args);

    // parameters of the non-bonded interactions
    args.check("cut", 1);
    args.check("epskap", 2);
    double cut, eps, kappa;
    {
      stringstream ss;
      ss << args["cut"];
      ss >> cut;
      ss.clear();
      ss.str("");
      Arguments::const_iterator start = args.lower_bound("epskap");
      Arguments::const_iterator end = args.upper_bound("epskap");
      for (; start != end; start++) {
        ss << start->second << endl;
      }
      ss >> eps >> kappa;
    }

    double cutmin = 0.0;
    if (args.count("cutmin") >= 1) {
      stringstream ss;
      ss << args["cutmin"];
      ss >> cutmin;
    }
    
    double cutrf = cut;
    if (args.count("cutrf") >= 1) {
      stringstream ss;
      ss << args["cutrf"];
      ss >> cutrf;
    }

    // create Time object to read the time from trajectory
    Time time(args);

    // add the atoms of groups A and B to the atom specifiers
    args.check("atomsA", 1);
    args.check("atomsB", 1);
    AtomSpecifier atomsA(sys);
    AtomSpecifier atomsB(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("atomsA");
      Arguments::const_iterator to = args.upper_bound("atomsA");
      for (; iter != to; iter++) {
        atomsA.addSpecifier(iter->second.c_str());
      }
      iter = args.lower_bound("atomsB");
      to = args.upper_bound("atomsB");
      for (; iter != to; iter++) {
        atomsB.addSpecifier(iter->second);
      }
    }
    // make sure the atoms are nicely sorted
    atomsA.sort();
    atomsB.sort();

    // get the parilist type
    args.check("pairlist", 1);
    string type = args["pairlist"];
    if (type != "ATOMIC" && type != "CHARGEGROUP") {
      throw gromos::Exception(argv[0], "ERROR: check argument of @pairlist: only ATOMIC or CHARGEGROUP is allowed");
    }
    // in case of pair-list type CHARGEGROUP, make sure that we have all complete charge groups in atomsB
    if (type == "CHARGEGROUP") {
      for (int b = 0; b < atomsB.size(); b++) {
        int molNum = atomsB.mol(b);
        int atomNum = atomsB.atom(b);
        int gromosNum = atomsB.gromosAtom(b);
        // get the atom and molecule number of the previous atom, if existing
        int prevMolNum = -1;
        int prevAtomNum = -1;
        if (molNum > 0 || atomNum > 0) {
          if (atomNum > 0) {
            prevMolNum = molNum;
            prevAtomNum = atomNum - 1;
          } else {
            prevMolNum = molNum - 1;
            prevAtomNum = sys.mol(prevMolNum).numAtoms() - 1;
          }
        }
        if (gromosNum > 0) {
          if (atomsB.sys()->mol(prevMolNum).topology().atom(prevAtomNum).chargeGroup() == 0) {
            stringstream ss;
            ss << "Pair-list algorithm CHARGEGROUP (@pairlist) together with incomplete charge groups within atomsB does not work. Check the charge group containing atom "
                    << gromosNum + 1 << " or use the pair-list algorithm ATOMIC (@pairlist).";
            throw gromos::Exception(argv[0], ss.str());
          }
        }
        while (atomsB.sys()->mol(molNum).topology().atom(atomNum).chargeGroup() == 0) {
          if (b + 1 >= atomsB.size()) {
            stringstream ss;
            ss << "Pair-list algorithm CHARGEGROUP (@pairlist) together with incomplete charge groups within atomsB does not work. Check the charge group containing atom "
                    << gromosNum + 1 << " or use the pair-list algorithm ATOMIC (@pairlist).";
            throw gromos::Exception(argv[0], ss.str());
          }
          if (gromosNum + 1 != atomsB.gromosAtom(b + 1)) {
            stringstream ss;
            ss << "Pair-list algorithm CHARGEGROUP (@pairlist) together with incomplete charge groups within atomsB does not work. Check the charge group containing atom "
                    << gromosNum + 1 << " or use the pair-list algorithm ATOMIC (@pairlist).";
            throw gromos::Exception(argv[0], ss.str());
          }
          b++;
          molNum = atomsB.mol(b);
          atomNum = atomsB.atom(b);
          gromosNum = atomsB.gromosAtom(b);
        }
        cerr << "Final Gromos Number: " << gromosNum << endl << endl;
      }
    }

    // in case of the 

    // print some information in case of @verbose
    if (args.count("verbose") >= 0) {
      cerr << "# cut = " << cut << endl;
      if (cutmin > 0.0) {
        cerr << "# cutmin = " << cutmin << endl;
      }
      if (args.count("cutrf") >= 1) {
        cerr << "# cutrf = " << cutrf << endl;
      }
      cerr << "# eps = " << eps << endl;
      cerr << "# kap = " << kappa << endl << "#" << endl;
      cerr << "# number of atoms in group A: " << atomsA.size() << endl;
      cerr << "# number of atoms in group B: " << atomsB.size() << endl << "#" << endl;
    }

    // loop over the trajectory files
    for (Arguments::const_iterator iter = args.lower_bound("trc");
            iter != args.upper_bound("trc"); ++iter) {

      // define input coordinates
      InG96 ic;
      ic.open(iter->second);
      ic.select("ALL");


      // initiate the pair list to be used later
      map<int, SimplePairlist> pl; // key (int) = gromos atom number
      //   for each atom of group A there is one key and one pair list
      // value (SimplePairlist): pair list from center
      //   atom "key" to all atoms within group B

      for (int a = 0; a < atomsA.size(); a++) {
        int gromosNumA = atomsA.gromosAtom(a);
        int molNumA = atomsA.mol(a);
        int atomNumA = atomsA.atom(a);
        SimplePairlist _pl(sys, *pbc, cut);
        pl[gromosNumA] = _pl;
        pl.find(gromosNumA)->second.setType(type);
        pl.find(gromosNumA)->second.setAtom(molNumA, atomNumA); // set the center atom
      }

      // a vector to keep the forces
      vector<Vec> force_CRF(atomsA.size());
      vector<Vec> force_LJ(atomsA.size());

      // the reaction field constant
      const double crf = ((2 - 2 * eps) * (1 + kappa * cut) - eps * (kappa * kappa * cut * cut)) /
              ((1 + 2 * eps)*(1 + kappa * cut) + eps * (kappa * kappa * cut * cut));

      // print the header of the table
      if (args.count("projvec") >= 1) {
        cout << "#" << setw(14) << "time" << setw(20) << "f_CRF_x" << setw(20) << "f_CRF_y" << setw(20) << "f_CRF_z" << setw(20)
                << "f_LJ_x" << setw(20) << "f_LJ_y" << setw(20) << "f_LJ_z" << setw(20)
                << "f_tot_x" << setw(20) << "f_tot_y" << setw(20) << "f_tot_z" << setw(20)
                << "projection_CRF" << setw(20) << "projection_LJ" << setw(20) << "projection_tot" << endl;
      } else {
        cout << "#" << setw(14) << "time" << setw(20) << "f_CRF_x" << setw(20) << "f_CRF_y" << setw(20) << "f_CRF_z" << setw(20)
                << "f_LJ_x" << setw(20) << "f_LJ_y" << setw(20) << "f_LJ_z" << setw(20)
                << "f_tot_x" << setw(20) << "f_tot_y" << setw(20) << "f_tot_z" << endl;
      }

      // loop over all frames
      while (!ic.eof()) {

        // read the configuration and the time of the current frame
        ic >> sys >> time;

        // the reference vector for the projection, if specified
        Vec e(0.0, 0.0, 0.0);
        if (args.count("projvec") >= 1) {
          VectorSpecifier vs(sys, pbc, args["projvec"]);
          e = vs().normalize();
        }

        // calculate the pair list for all atoms of group A
#ifdef OMP
#pragma omp parallel for
#endif
        for (int a = 0; a < atomsA.size(); a++) {

          // set the forces on atom a (within group A) to zero at the beginning
          force_CRF[a] = Vec(0.0, 0.0, 0.0);
          force_LJ[a] = Vec(0.0, 0.0, 0.0);

          int gromosNumA = atomsA.gromosAtom(a);
          Vec posA = atomsA.pos(a);
          double chargeA = atomsA.charge(a);
          map<int, SimplePairlist>::iterator it_pl = pl.find(gromosNumA); // the pair list to the current atom a within the group A
          it_pl->second.clear();
          it_pl->second.calc(atomsB, cutmin);
          it_pl->second.removeExclusions();

          // and the force acting on the atoms of group A by looping over
          // the according pair list

          for (int b = 0; b < it_pl->second.size(); b++) {
            //int molA = atomsA.mol(a);
            //int atomA = atomsA.atom(a);
            //int molB = it_pl->second.mol(b);
            //int atomB = it_pl->second.atom(b);
            //cerr << "# calculation interaction " << molA + 1 << ":" << atomA + 1 << "->" << molB + 1 << ":" << atomB + 1 << endl;
            Vec posB = pbc->nearestImage(posA, it_pl->second.pos(b), sys.box());
            Vec r_vec = posA - posB;
            double r = r_vec.abs();

            // the CRF force (without RF terms from excluded atoms)
            double chargeB = it_pl->second.charge(b);
            //cerr << "chargeA = " << chargeA << endl << "chargeB = " << chargeB << endl;
            double qq = chargeA * chargeB;
            Vec f = (qq * (physConst.get_four_pi_eps_i())) * (1 / (r * r * r) + crf / (cutrf * cutrf * cutrf)) * r_vec;
            force_CRF[a] += f;
            //cerr << "# qq = " << qq << endl;
            //cerr << "# f = (" << f[0] << "," << f[1] << "," << f[2] << ")\n";
            //cerr << "force_CRF[" << a << "] = (" << force_CRF[a][0] << "," << force_CRF[a][1] << "," << force_CRF[a][2] << ")\n";
            //cerr << "# crf = " << crf << endl;
            //cerr << "# cutrf = " << cutrf << endl;
            //cerr << "# r = " << r << endl << endl;

            // the LJ force
            LJType lj(gff.ljType(AtomPair(atomsA.iac(a), it_pl->second.iac(b))));
            // exclusions were already kicked out after building the pair list,
            // but we have to check for the special 1,4-interactions
            bool special14 = false;
            {
              // (14) exclusions are from the atom with the smaller gromos number to that one with the higher, so
              // make sure we agree with that convention ("smaller" atom excludes "bigger" atom, in terms of sequential numbers)
              int atomNum1 = atomsA.atom(a) <= it_pl->second.atom(b) ? atomsA.atom(a) : it_pl->second.atom(b);
              int atomNum2 = atomsA.atom(a) <= it_pl->second.atom(b) ? it_pl->second.atom(b) : atomsA.atom(a);
              // loop over the 14 exclusions of the "smaller" atom to see if the other one belong to the set of special 14 exclusion
              for (int e = 0; e < sys.mol(atomsA.mol(a)).topology().atom(atomNum1).exclusion14().size(); e++) {
                // if the two atoms are in a different molecule, there is nothing to do at all
                if (atomsA.mol(a) != it_pl->second.mol(b)) {
                  continue;
                }
                if (atomNum2 == sys.mol(atomsA.mol(a)).topology().atom(atomNum1).exclusion14().atom(e)) {
                  special14 = true;
                  break;
                }
              }
            }
            // get the LJ parameters, according to if it is a normal or special 14 LJ atom pair
            double c12, c6;
            if (special14) {
              c12 = lj.cs12();
              c6 = lj.cs6();
            } else {
              c12 = lj.c12();
              c6 = lj.c6();
            }
            //cerr << "c12 = " << c12 << endl << "c6 = " << c6 << endl;
            // and finally calculate the force of the LJ interaction
            double r2 = r*r;
            f = ((2 * c12 / (r2 * r2 * r2)) - c6) * (6 * r_vec / (r2 * r2 * r2 * r2));
            force_LJ[a] += f;
          }

          // add the RF terms from the excluded atoms to the force_CRF
          // (one may improve the code's speed, but at the moment we just want a working solution)
          for (int b = 0; b < atomsB.size(); b++) {
            // there is nothing to do if A and B are the same atom
            if (atomsA.gromosAtom(a) == atomsB.gromosAtom(b)) {
              continue;
            }
            // we have to loop over the exclusions of the atom with the lower
            // GROMOS number to see if the other one is in the set of excluded atoms
            bool excluded = false;
            if (atomsA.gromosAtom(a) < atomsB.gromosAtom(b)) {
              // loop over the exclusions of atom A to see if B is in this set
              int gromosNum1 = atomsA.gromosAtom(a);
              int gromosNum2 = atomsB.gromosAtom(b);
              for (int e = 0; e < atomsA.sys()->mol(atomsA.mol(a)).topology().atom(atomsA.atom(a)).exclusion().size(); e++) {
                int excludedGromosAtomNum = gromosNum1 + atomsA.sys()->mol(atomsA.mol(a)).topology().atom(atomsA.atom(a)).exclusion().atom(e) - atomsA.atom(a);
                if (excludedGromosAtomNum == gromosNum2) {
                  excluded = true;
                  break;
                }
              }
            } else {
              // loop over the exclusions of B to see if A is in this set
              int gromosNum1 = atomsB.gromosAtom(b);
              int gromosNum2 = atomsA.gromosAtom(a);
              for (int e = 0; e < atomsB.sys()->mol(atomsB.mol(b)).topology().atom(atomsB.atom(b)).exclusion().size(); e++) {
                int excludedGromosAtomNum = gromosNum1 + atomsB.sys()->mol(atomsB.mol(b)).topology().atom(atomsB.atom(b)).exclusion().atom(e) - atomsB.atom(b);
                if (excludedGromosAtomNum == gromosNum2) {
                  excluded = true;
                  break;
                }
              }
            }
            if (excluded) {
              double chargeB = atomsB.charge(b);
              double qq = chargeA * chargeB;
              Vec posB = pbc->nearestImage(posA, atomsB.pos(b), sys.box());
              Vec r_vec = posA - posB;
              Vec f = qq * (physConst.get_four_pi_eps_i()) * crf / (cutrf * cutrf * cutrf) * r_vec;
              force_CRF[a] += f;
            }
          }

        }

        // calculate the average CRF and LJ force over all atoms of group A
        Vec f_CRF(0.0, 0.0, 0.0);
        Vec f_LJ(0.0, 0.0, 0.0);
        for (int a = 0; a < atomsA.size(); a++) {
          f_CRF += force_CRF[a];
          f_LJ += force_LJ[a];
        }
        f_CRF /= (double) atomsA.size();
        f_LJ /= (double) atomsA.size();

        // in case there is a projection vector given , the output prints the force vector and its projection to e, otherwise
        cout.precision(9);
        if (e.abs2() > 0) {
          cout << setw(15) << time << scientific << setw(20)
                  << f_CRF[0] << setw(20) << f_CRF[1] << setw(20) << f_CRF[2] << setw(20)
                  << f_LJ[0] << setw(20) << f_LJ[1] << setw(20) << f_LJ[2] << setw(20)
                  << f_CRF[0] + f_LJ[0] << setw(20) << f_CRF[1] + f_LJ[1] << setw(20) << f_CRF[2] + f_LJ[2] << setw(20)
                  << f_CRF.dot(e) << setw(20) << f_LJ.dot(e) << setw(20) << f_CRF.dot(e) + f_LJ.dot(e) << endl;
        } else {
          cout << setw(15) << time << scientific << setw(20)
                  << f_CRF[0] << setw(20) << f_CRF[1] << setw(20) << f_CRF[2] << setw(20)
                  << f_LJ[0] << setw(20) << f_LJ[1] << setw(20) << f_LJ[2] << setw(20)
                  << f_CRF[0] + f_LJ[0] << setw(20) << f_CRF[1] + f_LJ[1] << setw(20) << f_CRF[2] + f_LJ[2] << endl;
        }

      }
    }

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

