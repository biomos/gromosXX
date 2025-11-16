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
 * @file shake_analysis.cc
 * Calculate all interactions for a selected set of atoms
 */

/**
 * @page programs Program Documentation
 *
 * @anchor shake_analysis
 * @section shake_analysis Calculate all interactions for a selected set of atoms
 * @author @ref co
 * @date 20.11.2004
 *
 * A SHAKE failure in one of the MD engines is one of the few indications that
 * something is going wrong in your simulation. Most often, there is a mismatch
 * between the topology and the set of coordinates, or an extremely high
 * interaction between particles is built up otherwise. Program shake_analysis
 * is a diagnostic tool that can be used to evaluate all interaction energies
 * for selected atoms, on a coordinate file right before or after a SHAKE
 * failure. The output can be limited by specifying the number of interactions
 * that should be displayed, or by giving an energy cutoff above which
 * interactions should be listed.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms" for which shake fails&gt; </td></tr>
 * <tr><td> \@cut</td><td>&lt;cut-off distance&gt; </td></tr>
 * <tr><td> [\@eps</td><td>&lt;epsilon for reaction field correction&gt;] </td></tr>
 * <tr><td> [\@kap</td><td>&lt;kappa for reaction field correction&gt;] </td></tr>
 * <tr><td> [\@top</td><td>&lt;number of non-bonded interactions per atom to print&gt;] </td></tr>
 * <tr><td> [\@higher</td><td>&lt;print energies higher than specified value&gt;] </td></tr>
 * <tr><td> [\@nocov</td><td>(do not print covalent interactions)] </td></tr>
 * <tr><td> \@coord</td><td>&lt;coordinate file&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  shake_analysis
    @topo    ex.top
    @pbc     r
    @atoms   1:34,35
    @coord   ex.coo
    @cut     1.4
    @eps     66.0
    @kap     0.0
    @top     10
    @higher  2e03
  @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <set>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "../src/gmath/Vec.h"
#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/gmath/Vec.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/BondType.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/AngleType.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/DihedralType.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/ImproperType.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/Property.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/SimplePairlist.h"
#include "../src/utils/Energy.h"
#include "../src/gromos/Exception.h"

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

using namespace std;

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "atoms" << "props" << "cut" << "eps" << "kap"
          << "top" << "coord" << "higher" << "nocov";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo    <molecular topology file>\n";
  usage += "\t@pbc     <boundary type>\n";
  usage += "\t@atoms   <atoms for which shake fails>\n";
  usage += "\t@cut     <cut-off distance>\n";
  usage += "\t[@eps    <epsilon for reaction field correction>]\n";
  usage += "\t[@kap    <kappa for reaction field correction>]\n";
  usage += "\t[@top    <number of non-bonded interactions per atom to print>]\n";
  usage += "\t[@higher <print energies higher than specified value>]\n";
  usage += "\t[@nocov  (do not print covalent interactions)]\n";
  usage += "\t@coord   <coordinate file>\n";


  try {
    Arguments args(argc, argv, knowns, usage);

    //  read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    GromosForceField gff(it.forceField());

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);

    // declare the energy class
    Energy en(sys, gff, *pbc);

    //  set atoms
    AtomSpecifier atoms(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("atoms");
      Arguments::const_iterator to = args.upper_bound("atoms");
      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        atoms.addSpecifier(spec);
      }
    }

    // get cutoff distance
    double cut = args.getValue<double>("cut", true);
    en.setCutOff(cut);


    // get RF variables
    double eps = args.getValue<double>("eps", false, 0.0);
    double kap = args.getValue<double>("kap", false, 0.0);
    en.setRF(eps, kap);

    // get top number
    int top = args.getValue<int>("top", false, -1);

    // get higher argument
    double higher = args.getValue<double>("higher", false, -1.0);

    // find the properties to calculate
    PropertyContainer props(sys, pbc);
    vector<int> num_prop(4, 0);
    vector<int> prop_types;

    for (int m = 0; m < sys.numMolecules(); m++) {
      BondIterator bi(sys.mol(m).topology());
      for (; bi; ++bi) {
        if (atoms.findAtom(m, bi()[0]) != -1 || atoms.findAtom(m, bi()[1]) != -1) {
          ostringstream os;
          os << "d%" << m + 1 << ":" << bi()[0] + 1 << "," << bi()[1] + 1;
          props.addSpecifier(os.str());
          num_prop[0]++;
          prop_types.push_back(bi().type());
        }
      }
    }
    for (int m = 0; m < sys.numMolecules(); m++) {
      AngleIterator ai(sys.mol(m).topology());
      for (; ai; ++ai) {
        if (atoms.findAtom(m, ai()[0]) != -1 ||
                atoms.findAtom(m, ai()[1]) != -1 ||
                atoms.findAtom(m, ai()[2]) != -1) {
          ostringstream os;
          os << "a%" << m + 1 << ":" << ai()[0] + 1 << "," << ai()[1] + 1 << ","
                  << ai()[2] + 1;
          props.addSpecifier(os.str());
          num_prop[1]++;
          prop_types.push_back(ai().type());
        }
      }
    }
    for (int m = 0; m < sys.numMolecules(); m++) {
      ImproperIterator ii(sys.mol(m).topology());
      for (; ii; ++ii) {
        if (atoms.findAtom(m, ii()[0]) != -1 ||
                atoms.findAtom(m, ii()[1]) != -1 ||
                atoms.findAtom(m, ii()[2]) != -1 ||
                atoms.findAtom(m, ii()[3]) != -1) {
          ostringstream os;
          os << "t%" << m + 1 << ":" << ii()[0] + 1 << "," << ii()[1] + 1 << ","
                  << ii()[2] + 1 << "," << ii()[3] + 1;
          props.addSpecifier(os.str());
          num_prop[2]++;
          prop_types.push_back(ii().type());
        }
      }
    }
    for (int m = 0; m < sys.numMolecules(); m++) {
      DihedralIterator di(sys.mol(m).topology());
      for (; di; ++di) {
        if (atoms.findAtom(m, di()[0]) != -1 ||
                atoms.findAtom(m, di()[1]) != -1 ||
                atoms.findAtom(m, di()[2]) != -1 ||
                atoms.findAtom(m, di()[3]) != -1) {
          ostringstream os;
          os << "t%" << m + 1 << ":" << di()[0] + 1 << "," << di()[1] + 1 << ","
                  << di()[2] + 1 << "," << di()[3] + 1;
          props.addSpecifier(os.str());
          num_prop[3]++;
          prop_types.push_back(di().type());

        }
      }
    }
    en.setProperties(props);

    // define input coordinate
    InG96 ic(args["coord"]);
    ic.select("ALL");
    ic >> sys;

    // we have to gather with any method to get covalent interactions
    // and charge-groups connected
    pbc->gathergr();

    // calculate the covalent energies
    if (args.count("nocov") < 0) {
      en.calcCov();

      // print out the covalent energies
      if (num_prop[0] + num_prop[1] + num_prop[2] + num_prop[3]) {
        cout << "--------------------------------------------------------------"
                << endl;

        cout << "Covalent interactions involving atoms ";
        for (unsigned int i = 0; i < atoms.size(); i++) {
          if (atoms.mol(i) < 0) cout << "s";
          else cout << atoms.mol(i) + 1;
          cout << ":" << atoms.atom(i) + 1 << " ";
        }
        cout << endl;
      }

      int count = 0;
      if (num_prop[0]) {
        cout << endl << "BONDS :" << endl << endl;
        cout << setw(4) << "mol"
                << setw(10) << "atom-"
                << setw(12) << "atom-"
                << setw(13) << "force-"
                << setw(10) << "b0"
                << setw(16) << "b in x"
                << setw(16) << "energy" << endl;
        cout << setw(4) << "# "
                << setw(10) << "numbers"
                << setw(12) << "names"
                << setw(13) << "constant"
                << endl;
      }
      for (int i = 0; i < num_prop[0]; i++, count++) {
        int type = prop_types[count];

        cout << setw(4) << props[count]->atoms().mol(0) + 1;
        cout << setw(5) << props[count]->atoms().atom(0) + 1 << "-"
                << setw(4) << props[count]->atoms().atom(1) + 1
                << setw(7) << sys.mol(props[count]->atoms().mol(0)).topology().
                atom(props[count]->atoms().atom(0)).name() << "-"
                << setw(4) << sys.mol(props[count]->atoms().mol(1)).topology().
                atom(props[count]->atoms().atom(1)).name();

        cout.precision(3);
        cout.setf(ios::scientific, ios::floatfield);

        cout << setw(13) << gff.bondType(type).fc();
        cout.setf(ios::fixed, ios::floatfield);
        cout << setw(10) << gff.bondType(type).b0();
        cout.precision(5);
        cout << setw(16) << props[count]->getValue();
        cout.setf(ios::scientific, ios::floatfield);
        cout << setw(16) << en.cov(count) << endl;
      }
      if (num_prop[1]) {
        cout << endl << "ANGLES :" << endl << endl;
        cout << setw(4) << "mol"
                << setw(15) << "atom-"
                << setw(17) << "atom-"
                << setw(13) << "force-"
                << setw(10) << "b0"
                << setw(16) << "b in x"
                << setw(16) << "energy" << endl;
        cout << setw(4) << "# "
                << setw(15) << "numbers"
                << setw(17) << "names"
                << setw(13) << "constant"
                << endl;
      }
      for (int i = 0; i < num_prop[1]; i++, count++) {
        int type = prop_types[count];
        cout << setw(4) << props[count]->atoms().mol(0) + 1;
        cout << setw(5) << props[count]->atoms().atom(0) + 1 << "-"
                << setw(4) << props[count]->atoms().atom(1) + 1 << "-"
                << setw(4) << props[count]->atoms().atom(2) + 1
                << setw(7) << sys.mol(props[count]->atoms().mol(0)).topology().
                atom(props[count]->atoms().atom(0)).name() << "-"
                << setw(4) << sys.mol(props[count]->atoms().mol(1)).topology().
                atom(props[count]->atoms().atom(1)).name() << "-"
                << setw(4) << sys.mol(props[count]->atoms().mol(2)).topology().
                atom(props[count]->atoms().atom(2)).name();
        cout.precision(3);
        cout.setf(ios::scientific, ios::floatfield);

        cout << setw(13) << gff.angleType(type).fc();
        cout.setf(ios::fixed, ios::floatfield);
        cout << setw(10) << gff.angleType(type).t0();
        cout.precision(5);
        cout << setw(16) << props[count]->getValue();
        cout.setf(ios::scientific, ios::floatfield);
        cout << setw(16) << en.cov(count) << endl;
      }
      if (num_prop[2]) {
        cout << endl << "IMPROPER DIHEDRALS :" << endl << endl;
        cout << setw(4) << "mol"
                << setw(20) << "atom-"
                << setw(22) << "atom-"
                << setw(13) << "force-"
                << setw(10) << "b0"
                << setw(16) << "b in x"
                << setw(16) << "energy" << endl;
        cout << setw(4) << "# "
                << setw(20) << "numbers"
                << setw(22) << "names"
                << setw(13) << "constant"
                << endl;
      }
      for (int i = 0; i < num_prop[2]; i++, count++) {
        int type = prop_types[count];

        cout << setw(4) << props[count]->atoms().mol(0) + 1;
        cout << setw(5) << props[count]->atoms().atom(0) + 1 << "-"
                << setw(4) << props[count]->atoms().atom(1) + 1 << "-"
                << setw(4) << props[count]->atoms().atom(2) + 1 << "-"
                << setw(4) << props[count]->atoms().atom(3) + 1
                << setw(7) << sys.mol(props[count]->atoms().mol(0)).topology().
                atom(props[count]->atoms().atom(0)).name() << "-"
                << setw(4) << sys.mol(props[count]->atoms().mol(1)).topology().
                atom(props[count]->atoms().atom(1)).name() << "-"
                << setw(4) << sys.mol(props[count]->atoms().mol(2)).topology().
                atom(props[count]->atoms().atom(2)).name() << "-"
                << setw(4) << sys.mol(props[count]->atoms().mol(3)).topology().
                atom(props[count]->atoms().atom(3)).name();
        cout.precision(3);
        cout.setf(ios::scientific, ios::floatfield);

        cout << setw(13) << gff.improperType(type).fc();
        cout.setf(ios::fixed, ios::floatfield);
        cout << setw(10) << gff.improperType(type).q0();
        cout.precision(5);
        cout << setw(16) << props[count]->getValue();
        cout.setf(ios::scientific, ios::floatfield);
        cout << setw(16) << en.cov(count) << endl;
      }
      if (num_prop[3]) {
        cout << endl << "DIHEDRAL ANGLES :" << endl << endl;
        cout << setw(4) << "mol"
                << setw(20) << "atom-"
                << setw(22) << "atom-"
                << setw(13) << "force-"
                << setw(6) << "pd"
                << setw(4) << "np"
                << setw(16) << "b in x"
                << setw(16) << "energy" << endl;
        cout << setw(4) << "# "
                << setw(20) << "numbers"
                << setw(22) << "names"
                << setw(13) << "constant"
                << endl;
      }
      for (int i = 0; i < num_prop[3]; i++, count++) {
        int type = prop_types[count];
        cout << setw(4) << props[count]->atoms().mol(0) + 1;
        cout << setw(5) << props[count]->atoms().atom(0) + 1 << "-"
                << setw(4) << props[count]->atoms().atom(1) + 1 << "-"
                << setw(4) << props[count]->atoms().atom(2) + 1 << "-"
                << setw(4) << props[count]->atoms().atom(3) + 1
                << setw(7) << sys.mol(props[count]->atoms().mol(0)).topology().
                atom(props[count]->atoms().atom(0)).name() << "-"
                << setw(4) << sys.mol(props[count]->atoms().mol(1)).topology().
                atom(props[count]->atoms().atom(1)).name() << "-"
                << setw(4) << sys.mol(props[count]->atoms().mol(2)).topology().
                atom(props[count]->atoms().atom(2)).name() << "-"
                << setw(4) << sys.mol(props[count]->atoms().mol(3)).topology().
                atom(props[count]->atoms().atom(3)).name();
        cout.precision(3);
        cout.setf(ios::scientific, ios::floatfield);

        cout << setw(13) << gff.dihedralType(type).fc();
        cout.setf(ios::fixed, ios::floatfield);
        cout.precision(1);
        cout << setw(6) << gff.dihedralType(type).pd();
        cout << setw(4) << gff.dihedralType(type).np();

        cout.precision(5);
        cout << setw(16) << props[count]->getValue();
        cout.setf(ios::scientific, ios::floatfield);
        cout << setw(16) << en.cov(count) << endl;
      }
    }

    double sum_e_vdw = 0.0, sum_e_crf = 0.0, sum_e_tot = 0.0;
    // double max_e_vdw = 0.0, max_e_crf = 0.0, max_e_tot = 0.0;

    // loop over the relevant atoms
    for (unsigned int i = 0; i < atoms.size(); i++) {

      //make a pairlist
      SimplePairlist pl(sys, *pbc, cut);
      pl.setAtom(atoms.mol(i), atoms.atom(i));
      pl.setType("CHARGEGROUP");
      pl.calc();
      pl.removeExclusions();

      int atom_i = pl.size();

      // add the atom itself to the pairlist
      pl.addAtom(atoms.mol(i), atoms.atom(i));

      // set these atoms for the energy class
      en.setAtoms(pl);

      // calculate all the interactions over the pairlist individually
      vector<double> vdw(atom_i), el(atom_i);
      for (int j = 0; j < atom_i; j++) {
        en.calcPair(j, atom_i, vdw[j], el[j]);
      }

      // filter out the largest absolute interactions
      set<int> inter;
      int max_inter, max_atom;
      if (top == -1 || top > atom_i) max_inter = atom_i;
      else max_inter = top;
      vector<int> order;
      double e_tot, max;

      for (int k = 0; k < max_inter; k++) {
        max = 0.0;
        max_atom = 0;

        for (int l = 0; l < atom_i; l++) {
          e_tot = fabs(vdw[l] + el[l]);
          if (e_tot >= max && inter.count(l) == 0) {
            max = e_tot;
            max_atom = l;
          }
        }
        if (higher > 0 && max < higher)
          break;

        inter.insert(max_atom);
        order.push_back(max_atom);
      }

      if (order.size() == 0) continue;

      // now write out the max_inter interactions
      cout << endl;
      cout << "--------------------------------------------------------------"
              << endl;
      if (max_inter != atom_i) cout << "Largest " << max_inter << " of ";
      cout << atom_i;
      cout << " non-bonded interactions with atom ";
      cout.precision(4);
      cout.setf(ios::right, ios::floatfield);
      if (pl.mol(atom_i) < 0) cout << "s";
      else cout << pl.mol(atom_i) + 1;
      cout << ":" << pl.atom(atom_i) + 1 << " (\"" << pl.name(atom_i)
              << "\"; IAC: " << pl.iac(atom_i) + 1 << "; charge: "
              << pl.charge(atom_i) << ")" << endl << endl;
      cout << setw(8) << "atom"
              << setw(5) << "name"
              << setw(5) << "IAC"
              << setw(11) << "charge"
              << setw(15) << "distance"
              << setw(15) << "vdw"
              << setw(15) << "coulomb"
              << setw(15) << "total" << endl;

      for (unsigned int k = 0; k < order.size(); k++) {
        cout.precision(4);
        cout.setf(ios::right, ios::floatfield);
        ostringstream os_k;
        if (pl.mol(order[k]) < 0) os_k << "s";
        else os_k << pl.mol(order[k]) + 1;
        os_k << ":" << pl.atom(order[k]) + 1;
        cout << setw(8) << os_k.str()
                << setw(5) << pl.name(order[k])
                << setw(5) << pl.iac(order[k]) + 1
                << setw(11) << pl.charge(order[k]);

        gmath::Vec v = pbc->nearestImage(*pl.coord(atom_i),
                *pl.coord(order[k]), sys.box());
        double d = (v - *pl.coord(atom_i)).abs();
        cout.precision(5);
        cout << setw(15) << d;
        cout.setf(ios::scientific, ios::floatfield);
        cout << setw(15) << vdw[order[k]]
                << setw(15) << el[order[k]]
                << setw(15) << vdw[order[k]] + el[order[k]];

        sum_e_vdw += vdw[order[k]];
        sum_e_crf += el[order[k]];
        sum_e_tot += vdw[order[k]] + el[order[k]];

        cout << endl;
      }
      cout << endl;

    }

    cout << "================================================================================\n"
            << setw(10) << " " << setw(20) << "vdw" << setw(20) << "crf" << setw(20) << "total" << "\n"
            << setw(10) << "total" << setw(20) << sum_e_vdw << setw(20) << sum_e_crf << setw(20) << sum_e_tot
            << "\n" << endl;

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}







