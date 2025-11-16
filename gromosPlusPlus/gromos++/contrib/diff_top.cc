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
 * @file diff_top.cc
 * lists parameter differences
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor diff_top
 * @section diff_top lists parameter differences
 * @author @ref ns
 * @date 07.07.2010
 *
 * Program diff_top can be used to list parameter differences between two topologies.
 * It will only consider the parameters which are actually used.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;the topology files&gt; </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
 diff_top
     @topo top1.top top2.top
   @endverbatim

 * <hr>
 */
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <map>
#include <string>

#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/System.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/BondType.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/AngleType.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/DihedralType.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/ImproperType.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/LJType.h"
#include "../src/gcore/LinearTopology.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/Constraint.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gromos/Exception.h"


using namespace args;
using namespace std;
using namespace gcore;
using namespace gio;
using namespace gmath;

static const double epsilon = 0.0000001;

int main(int argc, char** argv) {
  Argument_List knowns;
  knowns << "topo";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo        <topology files>\n";

  try {
    Arguments args(argc, argv, knowns, usage);
    if (args.count("topo") != 2) {
      throw gromos::Exception(argv[0], "This program takes two topologies!");
    }
    Arguments::const_iterator toparg = args.lower_bound("topo");
    InTopology it1(toparg->second);
    System sys1(it1.system());
    InTopology it2((++toparg)->second);
    System sys2(it2.system());
    GromosForceField ff1(it1.forceField());
    GromosForceField ff2(it2.forceField());

    if (fabs(ff1.fpepsi() - ff2.fpepsi()) > epsilon ||
        fabs(ff1.boltz() - ff2.boltz()) > epsilon ||
        fabs(ff1.hbar() - ff2.hbar()) > epsilon ||
        fabs(ff1.spdl() - ff2.spdl()) > epsilon) {
      cout << "Physical constants do not match: " << endl;
      cout.precision(8);
      cout << "FPEPSI   : " << setw(15) << ff1.fpepsi() << setw(15) << ff2.fpepsi() << endl;
      cout << "BOLTZ    : " << setw(15) << ff1.boltz() << setw(15) << ff2.boltz() << endl;
      cout << "HBAR     : " << setw(15) << ff1.hbar() << setw(15) << ff2.hbar() << endl;
      cout << "SPDL     : " << setw(15) << ff1.spdl() << setw(15) << ff2.spdl() << endl;
    }

    if (sys1.numMolecules() != sys2.numMolecules()) {
      cerr << "The number of molecules do not match. Cannot proceed." << endl;
      return 0;
    }

    map<int,int> iacs;
    for(int m = 0; m < sys1.numMolecules(); ++m) {
      const MoleculeTopology & mt1 = sys1.mol(m).topology();
      const MoleculeTopology & mt2 = sys2.mol(m).topology();
      if (mt1.numAtoms() != mt2.numAtoms()) {
        cerr << "The number of atoms in molecule " << m+1 << " do not match. Cannot proceed." << endl;
        return 0;
      }
      for(int a = 0; a < mt1.numAtoms(); ++a) {
        if (iacs.find(mt2.atom(a).iac()) == iacs.end()) {
          iacs[mt1.atom(a).iac()] = mt2.atom(a).iac();
        } else {
          if (iacs[mt1.atom(a).iac()] != mt2.atom(a).iac()) {
            cout << "Difference in atom " << a+1 << " of molecule " << m+1 << ":" << endl;
            cout << "IACs are not consistent: " << iacs[mt1.atom(a).iac()]+1 << " vs. " << mt2.atom(a).iac()+1 << endl;
          }
        }

        if (mt1.atom(a).name() != mt2.atom(a).name() ||
            mt1.atom(a).mass() != mt2.atom(a).mass() ||
            mt1.atom(a).chargeGroup() != mt2.atom(a).chargeGroup() ||
            mt1.atom(a).charge() != mt2.atom(a).charge()) {
          cout << "Difference in atom " << a+1 << " of molecule " << m+1 << ":" << endl;
          cout.precision(4);
          cout << setw(8) << mt1.atom(a).name() << setw(8) << mt1.atom(a).mass() << setw(8) << mt1.atom(a).charge() << setw(3) << mt1.atom(a).chargeGroup() << endl;
          cout << setw(8) << mt2.atom(a).name() << setw(8) << mt2.atom(a).mass() << setw(8) << mt2.atom(a).charge() << setw(3) << mt2.atom(a).chargeGroup() << endl;
        }
        bool exclusionsDiff = false;
        if (mt1.atom(a).exclusion().size() != mt2.atom(a).exclusion().size()) {
          exclusionsDiff = true;
        } else {
          for (int i = 0; i < mt1.atom(a).exclusion().size(); ++i) {
            if (!mt1.atom(a).exclusion().contains(mt2.atom(a).exclusion().atom(i))) {
              exclusionsDiff = true;
              break;
            }
          }
        }
        if (exclusionsDiff) {
          cout << "Difference in exclusions of atom " << a+1 << " of molecule " << m+1 << ":" << endl;
          for (int i = 0; i < mt1.atom(a).exclusion().size(); ++i) {
            cout << setw(5) << mt1.atom(a).exclusion().atom(i);
          }
          cout << endl;
          for (int i = 0; i < mt2.atom(a).exclusion().size(); ++i) {
            cout << setw(5) << mt2.atom(a).exclusion().atom(i);
          }
          cout << endl;
        }
        exclusionsDiff = false;
        if (mt1.atom(a).exclusion14().size() != mt2.atom(a).exclusion14().size()) {
          exclusionsDiff = true;
        } else {
          for (int i = 0; i < mt1.atom(a).exclusion14().size(); ++i) {
            int iaca1 = mt1.atom(a).iac(), iaca2 = mt2.atom(a).iac();
            if (mt1.atom(a).exclusion14().atom(i) != mt2.atom(a).exclusion14().atom(i)) {
              exclusionsDiff = true;
            } else {
              int iacb1 = mt1.atom(mt1.atom(a).exclusion14().atom(i)).iac();
              int iacb2 = mt2.atom(mt2.atom(a).exclusion14().atom(i)).iac();
              const LJType & t1 = ff1.ljType(AtomPair(iaca1,iacb1));
              const LJType & t2 = ff2.ljType(AtomPair(iaca2,iacb2));
              if (fabs(t1.cs6() - t2.cs6()) > epsilon ||
                fabs(t1.cs12() - t2.cs12()) > epsilon) {
                cout << "Difference in 1-4 pair " << m+1 << ":" << a+1 << "/" << m+1 << ":" << mt1.atom(a).exclusion14().atom(i)+1 << endl;
                cout.precision(5);
                cout << setw(14) << t1.cs6() << setw(14) << t1.cs12() << " vs. " << setw(14) << t2.cs6() << setw(14) << t2.cs12() << endl;
              }
            }
          }
        }
        if (exclusionsDiff) {
          cout << "Difference in 1,4-interactions of atom " << a+1 << " of molecule " << m+1 << ":" << endl;
          for (int i = 0; i < mt1.atom(a).exclusion14().size(); ++i) {
            cout << setw(5) << mt1.atom(a).exclusion14().atom(i);
          }
          cout << endl;
          for (int i = 0; i < mt2.atom(a).exclusion14().size(); ++i) {
            cout << setw(5) << mt2.atom(a).exclusion14().atom(i);
          }
          cout << endl;
        }
      }

      if (mt1.numBonds() == mt2.numBonds()) {
        BondIterator it1(mt1), it2(mt2);
        for(; it1 && it2; ++it1, ++it2) {
          if (it1()[0] != it2()[0] || 
              it1()[1] != it2()[1]) {
            cout << "Bond " << it1()[0] << "-" << it1()[1] << " vs. " << it2()[0] << "-" << it2()[1] << endl;
          } else {
            const BondType & b1 = ff1.bondType(it1().type());
            const BondType & b2 = ff2.bondType(it2().type());
            if (fabs(b1.b0() - b2.b0()) > epsilon ||
                fabs(b1.fc() - b2.fc()) > epsilon) {
              cout << "Bond Type " << b1.b0() << "/" << b1.fc() << " vs. " << b2.b0() << "/" << b2.fc() << endl;
            } 
          }
        }
      } else {
         cout << "The number of bonds does not match." << endl;
      }
      
      if (mt1.numAngles() == mt2.numAngles()) {
        AngleIterator it1(mt1), it2(mt2);
        for(; it1 && it2; ++it1, ++it2) {
          if (it1()[0] != it2()[0] || 
              it1()[1] != it2()[1] ||
              it1()[2] != it2()[2]) {
            cout << "Angle " << it1()[0] << "-" << it1()[1] << "-" << it1()[2] << " vs. " << it2()[0] << "-" << it2()[1] << "-" << it2()[2] << endl;
          } else {
            const AngleType & b1 = ff1.angleType(it1().type());
            const AngleType & b2 = ff2.angleType(it2().type());
            if (fabs(b1.t0() - b2.t0()) > epsilon ||
                fabs(b1.fc() - b2.fc()) > epsilon) {
              cout << "Angle Type " << b1.t0() << "/" << b1.fc() << " vs. " << b2.t0() << "/" << b2.fc() << endl;
            } 
          }
        }
      } else {
         cout << "The number of bond angles does not match." << endl;
      }

      if (mt1.numImpropers() == mt2.numImpropers()) {
        ImproperIterator it1(mt1), it2(mt2);
        for(; it1 && it2; ++it1, ++it2) {
          if (it1()[0] != it2()[0] ||
              it1()[1] != it2()[1] ||
              it1()[2] != it2()[2] ||
              it1()[3] != it2()[3]) {
            cout << "Improper "
                    << it1()[0] << "-" << it1()[1] << "-" << it1()[2] << "-" << it1()[3]
                    << " vs. "
                    << it2()[0] << "-" << it2()[1] << "-" << it2()[2] << "-" << it2()[3] << endl;
          } else {
            const ImproperType & b1 = ff1.improperType(it1().type());
            const ImproperType & b2 = ff2.improperType(it2().type());
            if (fabs(b1.q0() - b2.q0()) > epsilon ||
                fabs(b1.fc() - b2.fc()) > epsilon) {
              cout << "Improper Type " << b1.q0() << "/" << b1.fc() << " vs. " << b2.q0() << "/" << b2.fc() << endl;
            }
          }
        }
      } else {
         cout << "The number of improper dihedrals does not match." << endl;
      }

      if (mt1.numDihedrals() == mt2.numDihedrals()) {
        DihedralIterator it1(mt1), it2(mt2);
        for(; it1 && it2; ++it1, ++it2) {
          if (it1()[0] != it2()[0] ||
              it1()[1] != it2()[1] ||
              it1()[2] != it2()[2] ||
              it1()[3] != it2()[3]) {
            cout << "Dihedral "
                    << it1()[0] << "-" << it1()[1] << "-" << it1()[2] << "-" << it1()[3]
                    << " vs. "
                    << it2()[0] << "-" << it2()[1] << "-" << it2()[2] << "-" << it2()[3] << endl;
          } else {
            const DihedralType & b1 = ff1.dihedralType(it1().type());
            const DihedralType & b2 = ff2.dihedralType(it2().type());
            if (fabs(b1.fc() - b2.fc()) > epsilon ||
                fabs(b1.pdl() - b2.pdl()) > epsilon ||
                fabs(b1.np() - b2.np()) > epsilon) {
              cout << "Dihedral Type " 
                      << b1.fc() << "/" << b1.pdl() << "/" << b1.np()
                      << " vs. "
                      << b2.fc() << "/" << b2.pdl() << "/" << b2.np() << endl;
            }
          }
        }
      } else {
         cout << "The number of improper dihedrals does not match." << endl;
      }
    }

    // check atom pairs
    for(int m1 = 0; m1 < sys1.numMolecules(); ++m1) {
      for(int a1 = 0; a1 < sys1.mol(m1).topology().numAtoms(); ++a1) {
        for(int m2 = m1; m2 < sys1.numMolecules(); ++m2) {
          for(int a2 = a1; a2 < sys1.mol(m2).topology().numAtoms(); ++a2) {
            const LJType & t1 = ff1.ljType(AtomPair(sys1.mol(m1).topology().atom(a1).iac(),sys1.mol(m2).topology().atom(a2).iac()));
            const LJType & t2 = ff2.ljType(AtomPair(sys2.mol(m1).topology().atom(a1).iac(),sys2.mol(m2).topology().atom(a2).iac()));
            if (fabs(t1.c6() - t2.c6()) > epsilon || 
                fabs(t1.c12() - t2.c12()) > epsilon) {
              cout << "Difference in atompair " << m1+1 << ":" << a1+1 << "/" << m2+1 << ":" << a2+1 << endl;
              cout.precision(5);
              cout << setw(14) << t1.c6() << setw(14) << t1.c12() << " vs. " << setw(14) << t2.c6() << setw(14) << t2.c12() << endl;
            }
          }
        }
      }
    }

    // check solvent
    if (sys1.sol(0).topology().numAtoms() != sys2.sol(0).topology().numAtoms()) {
      cout << "Number of atoms in solvent topology does not match." << endl;
    } else {
      for(int a = 0; a < sys1.sol(0).topology().numAtoms(); ++a) {
        if (sys1.sol(0).topology().atom(a).name() != sys2.sol(0).topology().atom(a).name() ||
            sys1.sol(0).topology().atom(a).mass() != sys2.sol(0).topology().atom(a).mass() ||
            sys1.sol(0).topology().atom(a).charge() != sys2.sol(0).topology().atom(a).charge()) {
          cout << "Difference in atom " << a+1 << " of solvent."<< endl;
          cout.precision(4);
          cout << setw(8) << sys1.sol(0).topology().atom(a).name() << setw(8) << sys1.sol(0).topology().atom(a).mass() << setw(8) << sys1.sol(0).topology().atom(a).charge() << endl;
          cout << setw(8) << sys2.sol(0).topology().atom(a).name() << setw(8) << sys2.sol(0).topology().atom(a).mass() << setw(8) << sys2.sol(0).topology().atom(a).charge() << endl;
        }
      }
    }
    if (sys1.sol(0).topology().numConstraints() != sys2.sol(0).topology().numConstraints()) {
      cout << "Number of constraints in solvent topology does not match." << endl;
    } else {
      ConstraintIterator it1(sys1.sol(0).topology()), it2(sys1.sol(0).topology());
      for(; it1 && it2; ++it1, ++it2) {
        if(it1()[0] != it2()[0] ||
           it1()[1] != it2()[1] ||
           fabs(it1().dist() - it2().dist()) > epsilon) {
          cout << "Difference in solvent constraint:" << endl
                  << it1()[0] << "-" << it1()[1] << "/" << setw(12) << it1().dist() <<  endl
                  << it2()[0] << "-" << it2()[1] << "/" << setw(12) << it2().dist() <<  endl;
        }
      }
    }

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

