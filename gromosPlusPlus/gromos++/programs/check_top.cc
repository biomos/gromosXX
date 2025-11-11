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
 * @file check_top.cc
 * Check a topology for (consistency) errors
 */

/**
 * @page programs Program Documentation
 *
 * @anchor check_top
 * @section check_top Check a topology for (consistency) errors
 * @author @ref co
 * @date 7-6-07
 *
 * Making a correct molecular topology is one of the most important 
 * requirements for doing a successful simulation. check_top helps to remove 
 * often made errors from a topology in three ways. First, it runs some
 * standard tests on the molecular topology and warns if something unexpected
 * is observed in the topology. Second, it can calculate all bonded interaction
 * energies for a given set of coordinates, to determine the compatibility of
 * the topology with the coordinates. Third, it can check for consistency in
 * the force field parameters by comparing it to a specified set of building 
 * blocks and force field parameters.
 *
 * In the first phase check top tests that:
 * <ol>
 * <li> there is maximally one bond defined between any pair of atoms
 * <li> no atom appear
s twice in the definition of one given bond
 * <li> only bondtypes are used that are defined in the topology
 * <li> a bond angle is defined for the atoms involved in any two bonds sharing 
        one atom
 * <li> there is maximally one bond angle defined for a given set of three atoms
 * <li> atoms involved in a bond-angle definition are bound to the central atom
 * <li> no atom appears twice in the definition of one given bond angle
 * <li> only bond-angletypes are used that are defined in the topology
 * <li> an improper dihedral angle is defined centered on every atom that is
 *      boundto exactly three other atoms
 * <li> there is maximally one improper dihedral angle defined any set of four
        atoms
 * <li> atoms involved in an improper dihedral angle definition are bound
 * <li> no atom appears twice in the definition of one given improper dihedral
        angle
 * <li> only improper-dihedraltypes are used that are defined in the topology
 * <li> atoms involved in a proper dihedral angle are sequentially bound
 * <li> no atom appears twice in the definition of one given dihedral angle
 * <li> only dihedral-angletypes are used that are defined in the topology
 * <li> only atomtypes are used that are defined in the topology
 * <li> the sum of partial charges on atoms in one charge group is integer 
        valued
 * <li> excluded atoms are 1,2- or 1,3- or 1,4-neighbours
 * <li> atoms only have atoms with a higher sequence number in their exclusion
        list
 * <li> 1,2- or 1,3-neighbours are excluded
 * <li> 1,4-exclusions are separated by 3 bonds (1,4-neighbours)
 * <li> atoms only have atoms with a higher sequence number in their 
        1,4-exclusion list
 * <li> 1,4-neighbours are in the exclusion list or in the 1,4-exclusion list
 * <li> no exclusions or 1,4-exclusions are defined for the last atom in the
        topology
 * <li> the charge group code of the last atom in the topology is 1
 * </ol>
 *
 * Additionally, for atoms that are 1,4 neighbours but are listed as excluded
 * atoms a warning is written out that this is usually only the case if an
 * aromatic group is involved. Note that a topology that passes all these tests
 * is by no means guaranteed to be error-free. Conversely, some of these tests
 * are merely meant as warnings for the user which may point at errors in the
 * majority of cases. In some cases, the user may very well want to use a
 * topology that does not fulfill all tests.
 *
 * In the second phase, potential energies of all bonds, bond-angles, improper 
 * dihedral angles and proper dihedral angles are calculated and written out.
 * Abnormally large energies or deviations from ideal values may indicate an 
 * error in the topology, or an inconsistent set of coordinates with the
 * topology. See the program @ref shake_analysis for a similar check on the
 * non-bonded interaction energies.
 *
 * In the third phase check top can compare the topology with other building
 * blocks in a specified molecular topology building block file and the
 * corresponding interaction function parameter file. It checks if in the
 * molecular topology building block file we observe:
 * <ol>
 * <li> other atoms with the same name and the same integer atom code (IAC)
 * <li> other atoms with the specified IAC
 * <li> other atoms with the same IAC and mass
 * <li> other atoms with the same IAC and charge
 * <li> other bonds between atoms of the same IAC with the same bond type
 * <li> other bond angles between atoms of the same IAC with the same 
        bondangletype
 * <li> other improper dihedral angles between atoms of the same IAC with the
        same improper-dihedral type
 * <li> other dihedral angles between atoms of the same IAC with the same 
        dihedralangletype
 * </ol>
 *
 * In cases where the parameters specified in the program are not observed
 * anywhere else, or when they are not the most common parameter, the program
 * prints out a list of possible alternatives. Again, we stress that check_top
 * only points at possible inconsistencies and does not necessarily indicate
 * errors in your topology.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> [\@coord</td><td>&lt;coordinate file&gt; </td></tr>
 * <tr><td> [\@pbc</td><td>&lt;boundary type&gt; &lt;gather method&gt; </td></tr>
 * <tr><td> [\@build</td><td>&lt;building block file for consistency check&gt; </td></tr>
 * <tr><td> [\@param</td><td>&lt;parameter file for consistency check&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  check_top
    @topo      ex.top
    @coord     exref.coo
    @pbc       r
    @build     mtb45a3.dat
    @param     ifp45a3.dat
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/BondType.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/AngleType.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/ImproperType.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/CrossDihedral.h"
#include "../src/gcore/DihedralType.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/BuildingBlock.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InBuildingBlock.h"
#include "../src/gio/InParameter.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/Energy.h"
#include "../src/utils/CheckTopo.h"
#include "../src/utils/FfExpert.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "coord" << "build" << "param";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo     <topology>\n";
  usage += "\t[@coord   <coordinate file>]\n";
  usage += "\t[@pbc     <boundary type> <gather method>]\n";
  usage += "\t[@build   <building block file for consistency check>]\n";
  usage += "\t[@param   <parameter file for consistency check>]\n";



  try {
    Arguments args(argc, argv, knowns, usage);

    //  read topology
    InTopology it(args["topo"]);

    System sys(it.system());
    int nummol = sys.numMolecules();
    GromosForceField gff(it.forceField());

    System refSys(it.system());

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

    // declare the energy class
    Energy en(sys, gff, *pbc);

    bool do_energy = false;
    if (args.count("coord") > 0) do_energy = true;

    // set properties
    PropertyContainer props(sys, pbc);
    std::vector<int> numbonds(nummol);
    std::vector<int> numangles(nummol);
    std::vector<int> numimp(nummol);
    std::vector<int> numdih(nummol);
    std::vector<int> numcrossdih(nummol);
    int maxatomtype = 0, maxbondtype = 0, maxangletype = 0, maximptype = 0, maxdihtype = 0, maxcrossdihtype = 0;

    std::vector<double> totcharge(nummol);

    // loop over all bonds
    for (int m = 0; m < nummol; m++) {
      numbonds[m] = 0;
      BondIterator bi(sys.mol(m).topology());
      for (; bi; ++bi) {
        if (do_energy) {
          ostringstream os;
          os << "d%" << m + 1 << ":" << bi()[0] + 1 << "," << bi()[1] + 1;
          props.addSpecifier(os.str());
        }
        numbonds[m]++;
        if (bi().type() > maxbondtype) maxbondtype = bi().type();
      }
    }
    // loop over all angles
    for (int m = 0; m < nummol; m++) {
      numangles[m] = 0;
      AngleIterator ai(sys.mol(m).topology());
      for (; ai; ++ai) {
        if (do_energy) {
          ostringstream os;
          os << "a%" << m + 1 << ":" << ai()[0] + 1 << "," << ai()[1] + 1 << ","
                  << ai()[2] + 1;
          props.addSpecifier(os.str());
        }
        numangles[m]++;
        if (ai().type() > maxangletype) maxangletype = ai().type();
      }
    }
    // loop over all impropers
    for (int m = 0; m < nummol; m++) {
      numimp[m] = 0;
      ImproperIterator ii(sys.mol(m).topology());
      for (; ii; ++ii) {
        if (do_energy) {
          ostringstream os;
          os << "t%" << m + 1 << ":" << ii()[0] + 1 << "," << ii()[1] + 1 << ","
                  << ii()[2] + 1 << "," << ii()[3] + 1;
          props.addSpecifier(os.str());
        }
        numimp[m]++;
        if (ii().type() > maximptype) maximptype = ii().type();
      }
    }
    // loop over all dihedrals
    for (int m = 0; m < nummol; m++) {
      numdih[m] = 0;
      DihedralIterator di(sys.mol(m).topology());
      for (; di; ++di) {
        if (do_energy) {
          ostringstream os;
          os << "t%" << m + 1 << ":" << di()[0] + 1 << "," << di()[1] + 1 << ","
                  << di()[2] + 1 << "," << di()[3] + 1;
          props.addSpecifier(os.str());
        }
        numdih[m]++;
        if (di().type() > maxdihtype) maxdihtype = di().type();
      }
    }
    // loop over all cross dihedrals
    for (int m = 0; m < nummol; m++) {
      numcrossdih[m] = 0;
      CrossDihedralIterator di(sys.mol(m).topology());
      for (; di; ++di) {
        if (do_energy) {
          ostringstream os;
          os << "ct%" << m + 1 << ":" << di()[0] + 1 << "," << di()[1] + 1 << ","
                  << di()[2] + 1 << "," << di()[3] + 1 << "%" << m + 1 << ":"
                  << di()[4] + 1 << "," << di()[5] + 1 << ","
                  << di()[6] + 1 << "," << di()[7] + 1;
          props.addSpecifier(os.str());
        }
        numcrossdih[m]++;
        if (di().type() > maxcrossdihtype) maxcrossdihtype = di().type();
      }
    }

    // calculate the total charge
    for (int m = 0; m < nummol; m++) {
      totcharge[m] = 0.0;
      for (int a = 0; a < sys.mol(m).numAtoms(); a++) {
        totcharge[m] += sys.mol(m).topology().atom(a).charge();
        if (sys.mol(m).topology().atom(a).iac() > maxatomtype)
          maxatomtype = sys.mol(m).topology().atom(a).iac();
      }
    }

    if (args.count("coord") > 0) {

      // now, we are done preparing everything the real program starts here
      // calculate the values of all the properties
      // before reading coordinates????
      // props.calc();

      // parse them into the energy class
      en.setProperties(props);

      // define input coordinate
      InG96 ic;

      // open file
      ic.open(args["coord"]);
      ic.select("SOLUTE");

      // read coordinates and gather (gather method does not really matter)
      ic >> sys;
      (*pbc.*gathmethod)();

      // calculate the energies
      en.calc();
    }

    // That was it, now for the output
    int tnumbonds = 0, tnumangles = 0, tnumimp = 0, tnumdih = 0, tnumcrossdih = 0;
    double ttotcharge = 0.0;
    double chrg_precision = 0.0001;
    cout << "Topology contains " << sys.numMolecules() << " molecule";
    if (sys.numMolecules() > 1) cout << "s";
    cout << ":" << endl << endl;
    cout << setw(10) << "molecule"
            << setw(12) << "# atoms"
            << setw(12) << "# bonds"
            << setw(12) << "# angles"
            << setw(12) << "# impropers"
            << setw(12) << "# dihedrals"
            << setw(12) << "# cross-dihedrals"
            << setw(12) << "tot charge" << endl;
    for (int m = 0; m < nummol; m++) {
      if (fabs(totcharge[m]) < chrg_precision) totcharge[m] = 0.0;
      cout << setw(10) << m + 1
              << setw(12) << sys.mol(m).topology().numAtoms()
              << setw(12) << numbonds[m]
              << setw(12) << numangles[m]
              << setw(12) << numimp[m]
              << setw(12) << numdih[m]
              << setw(12) << numcrossdih[m]
              << setw(12) << totcharge[m] << endl;
      tnumbonds += numbonds[m];
      tnumangles += numangles[m];
      tnumimp += numimp[m];
      tnumdih += numdih[m];
      tnumcrossdih += numcrossdih[m];
      ttotcharge += totcharge[m];

    }

    // do some tests on the topology
    cout << endl
            << "Performing some basic checks on the charge groups, exclusions, "
            << "bonds, angles and improper dihedrals..." << endl;
    int error = 0;

    // loop over the molecules
    for (int m = 0; m < sys.numMolecules(); m++) {

      // define a topology check
      utils::CheckTopo ct(sys, m);
      ct.checkAll();
      error += ct.numErrors();
      if (ct.numErrors()) {
        cout << "--------------------" << endl;
        cout << "In Molecule " << m + 1 << ":" << endl << endl;
        for (int i = 0; i < ct.numErrors(); i++)
          cout << ct.error(i) << endl << endl;
        cout << "--------------------" << endl;
      }
    }
    // check whether all types are defined
    ostringstream par;
    int parerr = 0;
    if (maxatomtype >= gff.numAtomTypeNames()) {
      par << "Higher atom type used than defined: " << maxatomtype + 1
              << " > " << gff.numAtomTypeNames() << endl << endl;
      ++parerr;
    }
    if (maxbondtype >= gff.numBondTypes()) {
      par << "Higher bond type used than defined: " << maxbondtype + 1
              << " > " << gff.numBondTypes() << endl << endl;
      ++parerr;
    }
    if (maxangletype >= gff.numAngleTypes()) {
      par << "Higher angle type used than defined: " << maxangletype + 1
              << " > " << gff.numAngleTypes() << endl << endl;
      ++parerr;
    }
    if (maximptype >= gff.numImproperTypes()) {
      par << "Higher improper dihedral type used than defined: " << maximptype + 1
              << " > " << gff.numImproperTypes() << endl << endl;
      ++parerr;
    }
    if (maxdihtype >= gff.numDihedralTypes()) {
      par << "Higher dihedral type used than defined: " << maxdihtype + 1
              << " > " << gff.numDihedralTypes() << endl << endl;
      ++parerr;
    }
    if (maxcrossdihtype >= gff.numDihedralTypes()) {
      par << "Higher cross dihedral type used than defined: " << maxcrossdihtype + 1
              << " > " << gff.numDihedralTypes() << endl << endl;
      ++parerr;
    }
    error += parerr;
    if (parerr) {
      cout << "--------------------" << endl;
      cout << "In Parameters:" << endl << endl;
      cout << par.str();
      cout << "--------------------" << endl;
    }
    if (error == 0)
      cout << "ok" << endl;

    // possibly perform consistency check with the building block
    if (args.count("build") > 0) {
      if (args.count("param") <= 0)
        throw gromos::Exception("check_top",
              "For consistency check, give both @build and "
              "@param input flags");
      cout << "\n\nComparing parameters with other building blocks for "
              << "consistency\n"
              << "Building block file: " << args["build"]
              << "\nParameter file: " << args["param"] << endl;

      gio::InBuildingBlock ibb(args["build"]);
      gcore::BuildingBlock mtb(ibb.building());
      gio::InParameter ip(args["param"]);
      gcore::GromosForceField gffc(ip.forceField());

      utils::FfExpert ffexp(mtb);
      vector<utils::FfExpert::counter> v;

      // loop over all atoms
      for (int m = 0; m < nummol; ++m) {
        cout << "\n--------------------" << endl;
        cout << "In Molecule " << m + 1 << ":\n\n";

        cout << endl << sys.mol(m).numAtoms() << " ATOMS :" << endl << endl;
        cout << setw(8) << "atom-"
                << setw(8) << "atom-"
                << setw(6) << "IAC"
                << setw(9) << "mass"
                << setw(9) << "charge"
                << " charge group    consistency" << endl;
        cout << setw(8) << "number"
                << setw(8) << "name" << endl;

        for (int a = 0; a < sys.mol(m).numAtoms(); ++a) {
          bool inconsistency = false;
          ostringstream wrn;

          // to prepare, write the atom information to the warning
          wrn << setw(8) << a + 1 << ' '
                  << setw(7) << sys.mol(m).topology().atom(a).name()
                  << setw(6) << sys.mol(m).topology().atom(a).iac() + 1
                  << setw(9) << sys.mol(m).topology().atom(a).mass()
                  << setw(9) << sys.mol(m).topology().atom(a).charge()
                  << setw(4) << sys.mol(m).topology().atom(a).chargeGroup();


          // check the name with the IAC
          ffexp.name2iac(sys.mol(m).topology().atom(a).name().substr(0, 1), v);

          bool found = false;
          for (unsigned int i = 0; i < v.size(); ++i) {
            if (v[i].type == sys.mol(m).topology().atom(a).iac()) found = true;
          }
          if (!found) {
            if (!inconsistency)
              wrn << "             Inconsistency with building block file:\n";
            wrn << "\t\tNo atoms found with name "
                    << sys.mol(m).topology().atom(a).name() << " and IAC "
                    << sys.mol(m).topology().atom(a).iac() + 1
                    << "\n\t\tSuggested IACs (atomTypeName; occurence):\n";
            utils::sort(v, true);

            for (unsigned int i = 0; i < v.size(); ++i) {
              wrn << "\t\t\t"
                      << setw(4) << v[i].type + 1 << " ("
                      << setw(5) << gffc.atomTypeName(v[i].type) << "; "
                      << setw(4) << v[i].occurence << ")\n";
            }

            inconsistency = true;
          }

          // check the mass with the IAC
          ffexp.iac2mass(sys.mol(m).topology().atom(a).iac(), v);

          if (!v.size()) {
            if (!inconsistency)
              wrn << "             Inconsistency with building block file:\n";
            wrn << "\t\tNo atoms found with IAC "
                    << sys.mol(m).topology().atom(a).iac() + 1 << "\n"
                    << "\t\tMaximum IAC in parameter file: "
                    << gffc.numAtomTypeNames() << "\n";
            inconsistency = true;
          } else {

            found = false;
            for (unsigned int i = 0; i < v.size(); ++i) {
              if (gffc.findMass(v[i].type) == sys.mol(m).topology().atom(a).mass()) found = true;
            }
            if (!found) {
              if (!inconsistency)
                wrn << "             Inconsistency with building block file:\n";
              wrn << "\t\tNo atoms found with IAC "
                      << sys.mol(m).topology().atom(a).iac() + 1 << " and Mass "
                      << sys.mol(m).topology().atom(a).mass()
                      << "\n\t\tSuggested Masstype (mass; occurence):\n";
              utils::sort(v, true);

              for (unsigned int i = 0; i < v.size(); ++i) {
                wrn << "\t\t\t"
                        << setw(4) << v[i].type << " ("
                        << setw(9) << gffc.findMass(v[i].type) << "; " << setw(5)
                        << v[i].occurence << ")";
              }
              inconsistency = true;
            }

            // check the charge with the IAC
            ffexp.iac2charge(sys.mol(m).topology().atom(a).iac(), v);

            found = false;
            for (unsigned int i = 0; i < v.size(); ++i) {
              if (ffexp.charge(v[i].type) ==
                  sys.mol(m).topology().atom(a).charge()) found = true;
            }
            if (!found) {
              if (!inconsistency)
                wrn << "             Inconsistency with building block file:\n";
              wrn << "\t\tNo atoms found with IAC "
                      << sys.mol(m).topology().atom(a).iac() + 1 << " and charge "
                      << sys.mol(m).topology().atom(a).charge()
                      << "\n\t\tSuggested Charge (occurence):\n";

              utils::sort(v, false);
              for (unsigned int i = 0; i < v.size(); ++i) {
                wrn << "\t\t\t"
                        << setw(9) << ffexp.charge(v[i].type) << " ("
                        << setw(5) << v[i].occurence << ")\n";
              }
              inconsistency = true;
            }
          }
          if (!inconsistency) wrn << "             OK";
          wrn << endl;
          cout << wrn.str();

        }

        // now do the bonds
        cout << endl << numbonds[m] << " BONDS :" << endl << endl;
        cout << setw(4) << "mol"
                << setw(10) << "atom-"
                << setw(12) << "atom-"
                << setw(12) << "atom-"
                << setw(8) << "bond-"
                << setw(6) << "found"
                << setw(8) << "most"
                << "  alternatives" << endl;

        cout << setw(4) << "# "
                << setw(10) << "numbers"
                << setw(12) << "names"
                << setw(12) << "IAC"
                << setw(8) << "type"
                << setw(14) << "common"
                << endl;
        BondIterator bi(sys.mol(m).topology());

        for (; bi; ++bi) {
          int type = bi().type();
          cout << setw(4) << m + 1;
          cout << setw(5) << bi()[0] + 1 << "-" << setw(4) << bi()[1] + 1
                  << setw(7) << sys.mol(m).topology().atom(bi()[0]).name() << "-"
                  << setw(4) << sys.mol(m).topology().atom(bi()[1]).name()
                  << setw(7) << sys.mol(m).topology().atom(bi()[0]).iac() + 1 << "-"
                  << setw(4) << sys.mol(m).topology().atom(bi()[1]).iac() + 1
                  << setw(8) << type + 1;



          // create bond of the types
          Bond b(sys.mol(m).topology().atom(bi()[0]).iac(),
                 sys.mol(m).topology().atom(bi()[1]).iac(), 0);
          ffexp.iac2bond(b, v);

          bool found = false;
          bool first = false;
          utils::sort(v, false);

          for (unsigned int i = 0; i < v.size(); ++i) {
            if (v[i].type == type) {
              found = true;
              if (i == 0) first = true;
            }
          }

          int nalt = v.size();
          if (found) {
            cout << setw(6) << "X";
            --nalt;
          } else cout << setw(6) << " ";
          if (first) cout << setw(6) << "X";
          else cout << setw(6) << " ";
          cout << setw(6) << nalt << endl;
          if (nalt) {
            cout << "\t\t" << setw(4) << "type "
                    << setw(16) << "force constant"
                    << setw(14) << "bond length"
                    << setw(16) << "(occurrence)\n";

            for (unsigned int i = 0; i < v.size(); ++i) {
              cout << "\t\t"
                      << setw(4) << v[i].type + 1 << ": ";
              cout.precision(3);
              cout.setf(ios::scientific, ios::floatfield);

              cout << setw(15) << gffc.bondType(v[i].type).fc();
              cout.setf(ios::fixed, ios::floatfield);
              cout << setw(14) << gffc.bondType(v[i].type).b0();
              cout << "   ("
                      << setw(5) << v[i].occurence << ")\n";
            }
          }



        }
        // now do the angles
        cout << endl << numangles[m] << " ANGLES :" << endl << endl;
        cout << setw(4) << "mol"
                << setw(15) << "atom-"
                << setw(17) << "atom-"
                << setw(17) << "atom-"
                << setw(8) << "angle-"
                << setw(6) << "found"
                << setw(8) << "most"
                << "  alternatives" << endl;
        cout << setw(4) << "# "
                << setw(15) << "numbers"
                << setw(17) << "names"
                << setw(17) << "IAC"
                << setw(8) << "type"
                << setw(14) << "common"
                << endl;
        AngleIterator ai(sys.mol(m).topology());

        for (; ai; ++ai) {
          int type = ai().type();
          cout << setw(4) << m + 1;
          cout << setw(5) << ai()[0] + 1 << "-" << setw(4) << ai()[1] + 1
                  << "-" << setw(4) << ai()[2] + 1
                  << setw(7) << sys.mol(m).topology().atom(ai()[0]).name() << "-"
                  << setw(4) << sys.mol(m).topology().atom(ai()[1]).name() << "-"
                  << setw(4) << sys.mol(m).topology().atom(ai()[2]).name()
                  << setw(7) << sys.mol(m).topology().atom(ai()[0]).iac() + 1 << "-"
                  << setw(4) << sys.mol(m).topology().atom(ai()[1]).iac() + 1 << "-"
                  << setw(4) << sys.mol(m).topology().atom(ai()[2]).iac() + 1
                  << setw(8) << type + 1;

          // create angle of the types
          Angle a(sys.mol(m).topology().atom(ai()[0]).iac(),
                  sys.mol(m).topology().atom(ai()[1]).iac(),
                  sys.mol(m).topology().atom(ai()[2]).iac(), 0);
          ffexp.iac2angle(a, v);

          bool found = false;
          bool first = false;
          utils::sort(v, false);

          for (unsigned int i = 0; i < v.size(); ++i) {
            if (v[i].type == type) {
              found = true;
              if (i == 0) first = true;
            }
          }

          int nalt = v.size();
          if (found) {
            cout << setw(6) << "X";
            --nalt;
          } else cout << setw(6) << " ";
          if (first) cout << setw(6) << "X";
          else cout << setw(6) << " ";
          cout << setw(6) << nalt << endl;
          if (nalt) {
            cout << "\t\t" << setw(4) << "type "
                    << setw(16) << "force constant"
                    << setw(14) << "angle"
                    << setw(16) << "(occurrence)\n";

            for (unsigned int i = 0; i < v.size(); ++i) {
              cout << "\t\t"
                      << setw(4) << v[i].type + 1 << ": ";
              cout.precision(3);
              cout.setf(ios::scientific, ios::floatfield);

              cout << setw(15) << gffc.angleType(v[i].type).fc();
              cout.setf(ios::fixed, ios::floatfield);
              cout << setw(14) << gffc.angleType(v[i].type).t0();
              cout << "   ("
                      << setw(5) << v[i].occurence << ")\n";
            }
          }
        }

        // now do the impropers
        cout << endl << numimp[m] << " IMPROPER DIHEDRALS :" << endl << endl;
        cout << setw(4) << "mol"
                << setw(20) << "atom-"
                << setw(22) << "atom-"
                << setw(22) << "atom-"
                << setw(12) << "improper-"
                << setw(6) << "found"
                << setw(8) << "most"
                << "  alternatives" << endl;
        cout << setw(4) << "# "
                << setw(20) << "numbers"
                << setw(22) << "names"
                << setw(22) << "IAC"
                << setw(12) << "type"
                << setw(14) << "common"
                << endl;
        ImproperIterator ii(sys.mol(m).topology());

        for (; ii; ++ii) {
          int type = ii().type();
          cout << setw(4) << m + 1;
          cout << setw(5) << ii()[0] + 1 << "-" << setw(4) << ii()[1] + 1 << "-"
                  << setw(4) << ii()[2] + 1 << "-" << setw(4) << ii()[3] + 1
                  << setw(7) << sys.mol(m).topology().atom(ii()[0]).name() << "-"
                  << setw(4) << sys.mol(m).topology().atom(ii()[1]).name() << "-"
                  << setw(4) << sys.mol(m).topology().atom(ii()[2]).name() << "-"
                  << setw(4) << sys.mol(m).topology().atom(ii()[3]).name()
                  << setw(7) << sys.mol(m).topology().atom(ii()[0]).iac() + 1 << "-"
                  << setw(4) << sys.mol(m).topology().atom(ii()[1]).iac() + 1 << "-"
                  << setw(4) << sys.mol(m).topology().atom(ii()[2]).iac() + 1 << "-"
                  << setw(4) << sys.mol(m).topology().atom(ii()[3]).iac() + 1
                  << setw(12) << type + 1;

          // create improper of the types
          Improper i(sys.mol(m).topology().atom(ii()[0]).iac(),
                  sys.mol(m).topology().atom(ii()[1]).iac(),
                  sys.mol(m).topology().atom(ii()[2]).iac(),
                  sys.mol(m).topology().atom(ii()[3]).iac(), 0);
          ffexp.iac2improper(i, v);

          bool found = false;
          bool first = false;
          utils::sort(v, false);

          for (unsigned int i = 0; i < v.size(); ++i) {
            if (v[i].type == type) {
              found = true;
              if (i == 0) first = true;
            }
          }
          int nalt = v.size();
          if (found) {
            cout << setw(6) << "X";
            --nalt;
          } else cout << setw(6) << " ";
          if (first) cout << setw(6) << "X";
          else cout << setw(6) << " ";
          cout << setw(6) << nalt << endl;
          if (nalt) {
            cout << "\t\t" << setw(4) << "type "
                    << setw(16) << "force constant"
                    << setw(14) << "improper"
                    << setw(16) << "(occurrence)\n";

            for (unsigned int i = 0; i < v.size(); ++i) {
              cout << "\t\t"
                      << setw(4) << v[i].type + 1 << ": ";
              cout.precision(3);
              cout.setf(ios::scientific, ios::floatfield);

              cout << setw(15) << gffc.improperType(v[i].type).fc();
              cout.setf(ios::fixed, ios::floatfield);
              cout << setw(14) << gffc.improperType(v[i].type).q0();
              cout << "   ("
                      << setw(5) << v[i].occurence << ")\n";
            }
          }
        }
        // now do the dihedrals
        cout << endl << numdih[m] << " DIHEDRAL ANGLES :" << endl << endl;
        cout << setw(4) << "mol"
                << setw(20) << "atom-"
                << setw(22) << "atom-"
                << setw(22) << "atom-"
                << setw(12) << "dihedral-"
                << setw(6) << "found"
                << setw(8) << "most"
                << "  alternatives" << endl;
        cout << setw(4) << "# "
                << setw(20) << "numbers"
                << setw(22) << "names"
                << setw(22) << "IAC"
                << setw(12) << "type"
                << setw(14) << "common"
                << endl;
        DihedralIterator di(sys.mol(m).topology());

        for (; di; ++di) {
          int type = di().type();
          cout << setw(4) << m + 1;
          cout << setw(5) << di()[0] + 1 << "-" << setw(4) << di()[1] + 1 << "-"
                  << setw(4) << di()[2] + 1 << "-" << setw(4) << di()[3] + 1
                  << setw(7) << sys.mol(m).topology().atom(di()[0]).name() << "-"
                  << setw(4) << sys.mol(m).topology().atom(di()[1]).name() << "-"
                  << setw(4) << sys.mol(m).topology().atom(di()[2]).name() << "-"
                  << setw(4) << sys.mol(m).topology().atom(di()[3]).name()
                  << setw(7) << sys.mol(m).topology().atom(di()[0]).iac() + 1 << "-"
                  << setw(4) << sys.mol(m).topology().atom(di()[1]).iac() + 1 << "-"
                  << setw(4) << sys.mol(m).topology().atom(di()[2]).iac() + 1 << "-"
                  << setw(4) << sys.mol(m).topology().atom(di()[3]).iac() + 1
                  << setw(12) << type + 1;

          // create angle of the types
          Dihedral i(sys.mol(m).topology().atom(di()[0]).iac(),
                  sys.mol(m).topology().atom(di()[1]).iac(),
                  sys.mol(m).topology().atom(di()[2]).iac(),
                  sys.mol(m).topology().atom(di()[3]).iac(),
                  0);
          ffexp.iac2dihedral(i, v);

          bool found = false;
          bool first = false;
          utils::sort(v, false);

          for (unsigned int i = 0; i < v.size(); ++i) {
            if (v[i].type == type) {
              found = true;
              if (i == 0) first = true;
            }
          }
          int nalt = v.size();
          if (found) {
            cout << setw(6) << "X";
            --nalt;
          } else cout << setw(6) << " ";
          if (first) cout << setw(6) << "X";
          else cout << setw(6) << " ";
          cout << setw(6) << nalt << endl;
          if (nalt) {
            cout << "\t\t" << setw(4) << "type "
                    << setw(16) << "force constant"
                    << setw(14) << "phase shift"
                    << setw(14) << "multiplicity"
                    << setw(16) << "(occurrence)\n";

            for (unsigned int i = 0; i < v.size(); ++i) {
              cout << "\t\t"
                      << setw(4) << v[i].type + 1 << ": ";
              cout.precision(3);
              cout.setf(ios::scientific, ios::floatfield);

              cout << setw(15) << gffc.dihedralType(v[i].type).fc();
              cout.setf(ios::fixed, ios::floatfield);
              cout.precision(1);
              cout << setw(14) << gffc.dihedralType(v[i].type).pd();
              cout << setw(14) << gffc.dihedralType(v[i].type).np();

              cout << "   ("
                      << setw(5) << v[i].occurence << ")\n";
            }
          }
        }
      }

    }

    if (args.count("coord") > 0) {
      cout << endl << "Read in coordinates and calculated covalent energies:"
              << endl;

      int count = 0;
      int type = 0;
      std::vector<double> totbonds(nummol), totangles(nummol), totimp(nummol), totdih(nummol), totcrossdih(nummol);

      // loop over the properties once again to print
      // bonds
      cout << endl << tnumbonds << " BONDS :" << endl << endl;
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
      for (int m = 0; m < sys.numMolecules(); m++) {
        BondIterator bi(sys.mol(m).topology());
        totbonds[m] = 0;

        for (; bi; ++bi) {
          type = bi().type();
          cout << setw(4) << m + 1;
          cout << setw(5) << bi()[0] + 1 << "-" << setw(4) << bi()[1] + 1
                  << setw(7) << sys.mol(m).topology().atom(bi()[0]).name() << "-"
                  << setw(4) << sys.mol(m).topology().atom(bi()[1]).name();
          cout.precision(3);
          cout.setf(ios::scientific, ios::floatfield);

          cout << setw(13) << gff.bondType(type).fc();
          cout.setf(ios::fixed, ios::floatfield);
          cout << setw(10) << gff.bondType(type).b0();
          cout.precision(5);
          cout << setw(16) << props[count]->getValue();
          cout.setf(ios::scientific, ios::floatfield);
          cout << setw(16) << en.cov(count) << endl;
          totbonds[m] += en.cov(count);

          count++;

        }
      }
      // angles
      cout << endl << tnumangles << " ANGLES :" << endl << endl;
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
      for (int m = 0; m < sys.numMolecules(); m++) {
        totangles[m] = 0;
        AngleIterator bi(sys.mol(m).topology());

        for (; bi; ++bi) {
          type = bi().type();
          cout << setw(4) << m + 1;
          cout << setw(5) << bi()[0] + 1 << "-" << setw(4) << bi()[1] + 1
                  << "-" << setw(4) << bi()[2] + 1
                  << setw(7) << sys.mol(m).topology().atom(bi()[0]).name() << "-"
                  << setw(4) << sys.mol(m).topology().atom(bi()[1]).name() << "-"
                  << setw(4) << sys.mol(m).topology().atom(bi()[2]).name();
          cout.precision(3);
          cout.setf(ios::scientific, ios::floatfield);

          cout << setw(13) << gff.angleType(type).fc();
          cout.setf(ios::fixed, ios::floatfield);
          cout << setw(10) << gff.angleType(type).t0();
          cout.precision(5);
          cout << setw(16) << props[count]->getValue();
          cout.setf(ios::scientific, ios::floatfield);
          cout << setw(16) << en.cov(count) << endl;
          totangles[m] += en.cov(count);

          count++;

        }
      }
      // impropers
      cout << endl << tnumimp << " IMPROPER DIHEDRALS :" << endl << endl;
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
      for (int m = 0; m < sys.numMolecules(); m++) {
        totimp[m] = 0;
        ImproperIterator bi(sys.mol(m).topology());

        for (; bi; ++bi) {
          type = bi().type();
          cout << setw(4) << m + 1;
          cout << setw(5) << bi()[0] + 1 << "-" << setw(4) << bi()[1] + 1 << "-"
                  << setw(4) << bi()[2] + 1 << "-" << setw(4) << bi()[3] + 1
                  << setw(7) << sys.mol(m).topology().atom(bi()[0]).name() << "-"
                  << setw(4) << sys.mol(m).topology().atom(bi()[1]).name() << "-"
                  << setw(4) << sys.mol(m).topology().atom(bi()[2]).name() << "-"
                  << setw(4) << sys.mol(m).topology().atom(bi()[3]).name();
          cout.precision(3);
          cout.setf(ios::scientific, ios::floatfield);

          cout << setw(13) << gff.improperType(type).fc();
          cout.setf(ios::fixed, ios::floatfield);
          cout << setw(10) << gff.improperType(type).q0();
          cout.precision(5);
          cout << setw(16) << props[count]->getValue();
          cout.setf(ios::scientific, ios::floatfield);
          cout << setw(16) << en.cov(count) << endl;
          totimp[m] += en.cov(count);

          count++;

        }
      }
      // dihedrals
      cout << endl << tnumdih << " DIHEDRAL ANGLES :" << endl << endl;
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
      for (int m = 0; m < sys.numMolecules(); m++) {
        totdih[m] = 0;
        DihedralIterator bi(sys.mol(m).topology());
        vector<int> old(4, 0);
        for (; bi; ++bi) {
          type = bi().type();
          cout << setw(4) << m + 1;
          cout << setw(5) << bi()[0] + 1 << "-" << setw(4) << bi()[1] + 1 << "-"
                  << setw(4) << bi()[2] + 1 << "-" << setw(4) << bi()[3] + 1
                  << setw(7) << sys.mol(m).topology().atom(bi()[0]).name() << "-"
                  << setw(4) << sys.mol(m).topology().atom(bi()[1]).name() << "-"
                  << setw(4) << sys.mol(m).topology().atom(bi()[2]).name() << "-"
                  << setw(4) << sys.mol(m).topology().atom(bi()[3]).name();
          cout.precision(3);
          cout.setf(ios::scientific, ios::floatfield);

          cout << setw(13) << gff.dihedralType(type).fc();
          cout.setf(ios::fixed, ios::floatfield);
          cout.precision(1);
          cout << setw(6) << gff.dihedralType(type).pd();
          cout << setw(4) << gff.dihedralType(type).np();

          // don't give the energy if it is still the same dihedral.
          if (old[0] == bi()[0] && old[1] == bi()[1] && old[2] == bi()[2] && old[3] == bi()[3]) {
            cout << endl;
          } else {
            cout.precision(5);
            cout << setw(16) << props[count]->getValue();
            cout.setf(ios::scientific, ios::floatfield);
            cout << setw(16) << en.cov(count) << endl;
            totdih[m] += en.cov(count);
          }

          // set the old dihedral
          for (unsigned int i = 0; i < 4; ++i) old[i] = bi()[i];

          count++;

        }
      }
      // cross dihedrals
      cout << endl << tnumcrossdih << " CROSS-DIHEDRAL ANGLES :" << endl << endl;
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
      for (int m = 0; m < sys.numMolecules(); m++) {
        totcrossdih[m] = 0;
        CrossDihedralIterator bi(sys.mol(m).topology());
        vector<int> old(8, 0);
        for (; bi; ++bi) {
          type = bi().type();
          cout << setw(4) << m + 1;
          cout << setw(5) << bi()[0] + 1 << "-" << setw(4) << bi()[1] + 1 << "-"
                  << setw(4) << bi()[2] + 1 << "-" << setw(4) << bi()[3] + 1
                  << setw(5) << bi()[4] + 1 << "-" << setw(4) << bi()[5] + 1 << "-"
                  << setw(4) << bi()[6] + 1 << "-" << setw(4) << bi()[7] + 1
                  << setw(7) << sys.mol(m).topology().atom(bi()[0]).name() << "-"
                  << setw(4) << sys.mol(m).topology().atom(bi()[1]).name() << "-"
                  << setw(4) << sys.mol(m).topology().atom(bi()[2]).name() << "-"
                  << setw(4) << sys.mol(m).topology().atom(bi()[3]).name() << "-"
                  << setw(7) << sys.mol(m).topology().atom(bi()[4]).name() << "-"
                  << setw(4) << sys.mol(m).topology().atom(bi()[5]).name() << "-"
                  << setw(4) << sys.mol(m).topology().atom(bi()[6]).name() << "-"
                  << setw(4) << sys.mol(m).topology().atom(bi()[7]).name();
          cout.precision(3);
          cout.setf(ios::scientific, ios::floatfield);

          cout << setw(13) << gff.dihedralType(type).fc();
          cout.setf(ios::fixed, ios::floatfield);
          cout.precision(1);
          cout << setw(6) << gff.dihedralType(type).pd();
          cout << setw(4) << gff.dihedralType(type).np();

          // don't give the energy if it is still the same dihedral.
          if (old[0] == bi()[0] && old[1] == bi()[1] && old[2] == bi()[2] && old[3] == bi()[3] &&
              old[4] == bi()[4] && old[5] == bi()[5] && old[6] == bi()[6] && old[7] == bi()[7] ) {
            cout << endl;
          } else {
            cout.precision(5);
            cout << setw(16) << props[count]->getValue();
            cout.setf(ios::scientific, ios::floatfield);
            cout << setw(16) << en.cov(count) << endl;
            totcrossdih[m] += en.cov(count);
          }

          // set the old dihedral
          for (unsigned int i = 0; i < 8; ++i) old[i] = bi()[i];

          count++;

        }
      }
      // now summarize some energies
      cout.setf(ios::scientific, ios::floatfield);
      double ttbonds = 0, ttangles = 0, ttimp = 0, ttdih = 0, ttcrossdih = 0;

      cout << endl << "SUMMARY :" << endl << endl;
      cout << "Total energies" << endl;
      cout << setw(10) << "molecule"
              << setw(15) << "bonds"
              << setw(15) << "angles"
              << setw(15) << "impropers"
              << setw(15) << "dihedrals"
              << setw(15) << "crossdih."
              << setw(15) << "total" << endl;
      for (int m = 0; m < nummol; m++) {
        cout << setw(10) << m + 1
                << setw(15) << totbonds[m]
                << setw(15) << totangles[m]
                << setw(15) << totimp[m]
                << setw(15) << totdih[m]
                << setw(15) << totcrossdih[m]
                << setw(15) << totbonds[m] + totangles[m] + totimp[m] + totdih[m] + totcrossdih[m] << endl;
        ttbonds += totbonds[m];
        ttangles += totangles[m];
        ttimp += totimp[m];
        ttdih += totdih[m];
        ttcrossdih += totcrossdih[m];
      }
      cout << endl;

      if (nummol > 1)
        cout << setw(10) << "total"
        << setw(15) << ttbonds
              << setw(15) << ttangles
              << setw(15) << ttimp
              << setw(15) << ttdih
              << setw(15) << ttcrossdih
              << setw(15) << ttbonds + ttangles + ttimp + ttdih + ttcrossdih << endl;
      cout << setw(10) << "average"
              << setw(15) << ttbonds / tnumbonds
              << setw(15) << ttangles / tnumangles
              << setw(15) << ttimp / tnumimp
              << setw(15) << ttdih / tnumdih
              << setw(15) << ttcrossdih / tnumcrossdih<< endl;

      const double totcov = ttbonds + ttangles + ttimp + ttdih + ttcrossdih;

      cout << endl << "Total covalent energy: " << totcov << endl;

      if (error)
        cout << endl << "There were " << error << " warnings about the charge"
              << " groups, exclusions, bonds, angles or improper dihedrals"
              << endl;

    }
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}







