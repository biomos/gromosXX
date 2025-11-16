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

// gio_OutTopology.cc
#include "OutTopology.h"

#include <cassert>
#include <set>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "../gcore/BondType.h"
#include "../gcore/Bond.h"
#include "../gcore/AngleType.h"
#include "../gcore/Angle.h"
#include "../gcore/DihedralType.h"
#include "../gcore/Dihedral.h"
#include "../gcore/CrossDihedral.h"
#include "../gcore/ImproperType.h"
#include "../gcore/Improper.h"
#include "../gcore/LJType.h"
#include "../gcore/LJException.h"
#include "../gcore/LJExceptionType.h"
#include "../gcore/CGType.h"
#include "../gcore/AtomPair.h"
#include "../gcore/Exclusion.h"
#include "../gcore/Constraint.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Molecule.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/System.h"
#include "../gcore/GromosForceField.h"
#include "../gcore/VirtualAtomType.h"
#include "../gcore/VirtualAtoms.h"
#include "../args/Arguments.h"
#include "../utils/AtomSpecifier.h"
#include "../utils/VirtualAtom.h"
#include "../gromos/Exception.h"

using gio::OutTopology;
using namespace gcore;
using namespace std;

OutTopology::OutTopology(ostream &os) : d_os(os) {
}

OutTopology::~OutTopology() {
}

void OutTopology::setTitle(const string &title) {
  d_title = title;
}

//OutTopology &OutTopology::operator<<(const gcore::Simulation &sim){

void OutTopology::write(const gcore::System &sys, const gcore::GromosForceField &gff) {

  // Title block
  d_os << "TITLE\n" << d_title << "\nEND\n";


  if (args::Arguments::outG96) {
    // TOPPHYSCON block
    d_os.precision(10);

    d_os << "TOPPHYSCON\n"
            << "# FPEPSI: 1.0/(4.0*PI*EPS0) (EPS0 is the permittivity of vacuum)\n"
            << gff.fpepsi()
            << "\n# HBAR: Planck's constant HBAR = H/(2* PI)\n"
            << gff.hbar()
            << "\nEND\n";
  } else {
    // PHYSICALCONSTANTS block
    d_os.precision(10);

    d_os << "PHYSICALCONSTANTS\n"
            << "# FPEPSI: 1.0/(4.0*PI*EPS0) (EPS0 is the permittivity of vacuum)\n"
            << gff.fpepsi()
            << "\n# HBAR: Planck's constant HBAR = H/(2* PI)\n"
            << gff.hbar()
            << "\n# SPDL: Speed of light (nm/ps)\n"
            << gff.spdl()
            << "\n# BOLTZ: Boltzmann's constant kB\n"
            << gff.boltz()
            << "\nEND\n";
  }

  // TOPVERSION block
  if (args::Arguments::outG96) {
    d_os << "TOPVERSION\n1.7\nEND\n";
  } else {
    d_os << "TOPVERSION\n2.0\nEND\n";
  }

  // ATOMTYPENAME block
  d_os << "ATOMTYPENAME\n"
          << "# NRATT: number of van der Waals atom types\n";
  int num = gff.numAtomTypeNames();
  d_os << num << "\n";
  d_os << "# TYPE: atom type names\n";
  for (int i = 0; i < num; ++i) {
    if (i > 0 && !(i % 10))d_os << "# " << i << "\n";
    d_os << gff.atomTypeName(i) << "\n";
  }
  d_os << "END\n";

  // RESNAME block
  d_os << "RESNAME\n"
          << "# NRAA2: number of residues in a solute molecule\n";
  num = 0;
  for (int i = 0; i < sys.numMolecules(); ++i)
    num += sys.mol(i).topology().numRes();

  d_os << num << "\n"
          << "# AANM: residue names\n";
  for (int i = 0, count = 0; i < sys.numMolecules(); ++i)
    for (int j = 0; j < sys.mol(i).topology().numRes(); ++j, ++count) {
      if (count > 0 && !(count % 10))d_os << "# " << count << "\n";
      d_os << sys.mol(i).topology().resName(j) << "\n";
    }
  d_os << "END\n";

  // SOLUTEATOM block
  d_os << "SOLUTEATOM\n"
          << "#   NRP: number of solute atoms\n";
  num = 0;
  for (int i = 0; i < sys.numMolecules(); ++i)
    num += sys.mol(i).numAtoms();

  d_os << setw(5) << num << "\n";
  d_os << "#  ATNM: atom number\n"
          << "#  MRES: residue number\n"
          << "#  PANM: atom name of solute atom\n"
          << "#   IAC: integer (van der Waals) atom type code\n"
          << "#  MASS: mass of solute atom\n"
          << "#    CG: charge of solute atom\n"
          << "#   CGC: charge group code (0 or 1)\n"
          << "#   INE: number of excluded atoms\n"
          << "# INE14: number of 1-4 interactions\n"
          << "# ATNM MRES PANM IAC     MASS       CG  CGC INE\n"
          << "#                                           INE14\n";

  for (int i = 0, offatom = 1, offres = 1; i < sys.numMolecules(); ++i) {
    for (int j = 0; j < sys.mol(i).numAtoms(); ++j) {
      d_os.precision(5);
      d_os.setf(ios::fixed, ios::floatfield);
      d_os << setw(6) << offatom + j 
           << ' ' << setw(4) << sys.mol(i).topology().resNum(j) + offres
           << ' ' << setw(4) << sys.mol(i).topology().atom(j).name()
           << ' ' << setw(3) << sys.mol(i).topology().atom(j).iac() + 1
           << ' ' << setw(8) << sys.mol(i).topology().atom(j).mass()
           << ' ' << setw(8) << sys.mol(i).topology().atom(j).charge()
           << ' ' << setw(2) << sys.mol(i).topology().atom(j).chargeGroup();
      // Exclusions
      d_os << ' ' << setw(5) << sys.mol(i).topology().atom(j).exclusion().size();
      for (int k = 0; k < sys.mol(i).topology().atom(j).exclusion().size(); ++k) {
        if (k % 6 == 0 && k != 0)
          d_os << "\n"
                << "                                               ";
        d_os << ' ' << setw(5) << sys.mol(i).topology().atom(j).exclusion().atom(k) + offatom;
      }
      d_os << "\n"
              << "                                              "
              << sys.mol(i).topology().atom(j).exclusion14().size();
      for (int k = 0; k < sys.mol(i).topology().atom(j).exclusion14().size(); ++k) {
        if (k % 6 == 0 && k != 0)
          d_os << "\n"
                << "                                             ";
        d_os << ' ' << setw(5) << sys.mol(i).topology().atom(j).exclusion14().atom(k) + offatom;
      }

      d_os << "\n";
    }
    offres += sys.mol(i).topology().numRes();
    offatom += sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  // SOLUTEPOLARISATION block
  // check whether we have polarisable solute atoms;
  num = 0;
  for (int i = 0; i < sys.numMolecules(); ++i) {
    for (int j = 0; j < sys.mol(i).numAtoms(); ++j) {
      if (sys.mol(i).topology().atom(j).isPolarisable()) ++num;
    }
  }
  if (num) {
    d_os << "SOLUTEPOLARISATION\n"
            << "# NPPOL: number of polarisable solute atoms\n"
            << setw(5) << num << "\n"
            << "# IPOLP: atom sequence number of the polarisable solute atom\n"
            << "#  ALPP: polarisability of the solute atom\n"
            << "# QPOLP: size of charge-on-spring connected to polarisable solute atoms\n"
            << "# ENOTP: damping level for polarisation\n"
            << "#   EPP: damping power for polarisation\n"
            << "#  GAMP: gamma of for the off-site pol. centre construction\n"
            << "#    IP: first atom for the off-site pol. centre construction\n"
            << "#    JP: second atom for the off-site pol. centre construction\n#\n"
            << "# IPOLP" << setw(14) << "ALPP" << setw(15) << "QPOLP" << setw(15)
            << "ENOTP" << setw(15) << "EPP" << setw(15) << "GAMP"
            << setw(15) << "IP" << setw(15) << "JP" << endl;
    for (int i = 0, offatom = 1; i < sys.numMolecules(); offatom += sys.mol(i++).numAtoms()) {
      for (int j = 0; j < sys.mol(i).numAtoms(); ++j) {
        if (!sys.mol(i).topology().atom(j).isPolarisable()) continue;
//        std::cout <<" sys.mol(i).topology().atom(j).poloffsiteI(): " << sys.mol(i).topology().atom(j).poloffsiteI() << "\n sys.mol(i).topology().atom(j).poloffsiteJ(): " << sys.mol(i).topology().atom(j).poloffsiteJ() << "\n offatom: " << offatom << std::endl;
        d_os.precision(7);
        d_os.setf(ios::fixed, ios::floatfield);
        d_os << setw(5) << offatom + j << ' '
                << setw(15) << sys.mol(i).topology().atom(j).polarisability()
                << setw(15) << sys.mol(i).topology().atom(j).cosCharge()
                << setw(15) << sys.mol(i).topology().atom(j).dampingLevel()
                << setw(15) << sys.mol(i).topology().atom(j).dampingPower()
                << setw(15) << sys.mol(i).topology().atom(j).poloffsiteGamma()
                << setw(15) << sys.mol(i).topology().atom(j).poloffsiteI() + offatom 
                << setw(15) << sys.mol(i).topology().atom(j).poloffsiteJ() + offatom
                << "\n";
      }
    }
    d_os << "END\n";
  }

  // CGSOLUTE block
  // check whether we have coarse grained solute atoms;
  num = 0;
  for (int i = 0; i < sys.numMolecules(); ++i) {
    for (int j = 0; j < sys.mol(i).numAtoms(); ++j) {
      if (sys.mol(i).topology().atom(j).isCoarseGrained()) ++num;
    }
  }
  if (num) {
    d_os << "CGSOLUTE\n"
            << "# NCGB: number of coarse grained particle ranges\n";
    unsigned int tot_atoms = 0;
    for (int i = 0, offatom = 1; i < sys.numMolecules(); offatom += sys.mol(i++).numAtoms()) {
      for (int j = 0; j < sys.mol(i).numAtoms(); ++j) {
        ++tot_atoms;
      }
    }
    vector<vector<unsigned int> > ranges;
    vector<unsigned int> range(3);
    // get the coarse grained range(s) of particles
    range[0] = tot_atoms;
    for (int i = 0, offatom = 1; i < sys.numMolecules(); offatom += sys.mol(i++).numAtoms()) {
      for (int j = 0; j < sys.mol(i).numAtoms(); ++j) {
        if (sys.mol(i).topology().atom(j).isCoarseGrained()) {
          if ((offatom + j) < int(range[0])) range[0] = offatom + j;
          if (j == (sys.mol(i).numAtoms() - 1)) { // last atom?
            if (i < (sys.numMolecules() - 1)) { // not the last molecule
              if (!sys.mol(i + 1).topology().atom(0).isCoarseGrained() // next molecule not CG?
                  || sys.mol(i).topology().atom(j).cg_factor() !=
                  sys.mol(i + 1).topology().atom(0).cg_factor()) { // or different CG type?
                range[1] = offatom + j;
                range[2] = sys.mol(i).topology().atom(j).cg_factor();
                ranges.push_back(range);
                range[0] = tot_atoms;
              }
            } else { // of the last molecules
              range[1] = offatom + j;
              range[2] = sys.mol(i).topology().atom(j).cg_factor();
              ranges.push_back(range);
            }
          } else {
            if (!sys.mol(i).topology().atom(j + 1).isCoarseGrained() // next particle in molecule not CG?
                || sys.mol(i).topology().atom(j).cg_factor() !=
                  sys.mol(i).topology().atom(j + 1).cg_factor()) { // or different CG type?
              range[1] = offatom + j;
              range[2] = sys.mol(i).topology().atom(j).cg_factor();
              ranges.push_back(range);
              range[0] = tot_atoms;
            }
          }
        } // if j coarse grained
      } // for loop atoms in molecule
    } // for loop molecules
    d_os << setw(5) << ranges.size() << "\n"
            << "# NRCGF[1..NCGB]: sequence number of the first coarse grained solute particle in range\n"
            << "# NRCGL[1..NCGB]: sequence number of the last coarse grained solute particle in range\n"
            << "# MSCAL[1..NCGB]: scaling factor for pressure correction\n"
            << "# NRCGF     NRCGL     MSCAL\n";
    for (unsigned int i = 0; i < ranges.size(); ++i) {
      d_os.precision(7);
      d_os.setf(ios::fixed, ios::floatfield);
      d_os << setw(5) << ranges[i][0] << "   " << setw(5) << ranges[i][1]
           << setw(10) << ranges[i][2] << "\n";
    }

    d_os << "END\n";
  }

  if(gff.numVirtualAtomTypes()){
    d_os << "VIRTUALATOMTYPE\n"
	 << "# NVAT\n"
         << setw(6) << gff.numVirtualAtomTypes() << "\n"
         << "#  VAT   DIS1   DIS2\n";
    int num = gff.numVirtualAtomTypes();
    for(int i =0; i< num; i++){
      d_os << setw(6) << gff.virtualAtomTypeLine(i).code()
           << setw(9) << gff.virtualAtomTypeLine(i).dis1()
           << setw(9) << gff.virtualAtomTypeLine(i).dis2() << endl;
    }
    d_os << "END\n";

  }
  if(sys.vas().numVirtualAtoms()){
    d_os << "VIRTUALATOMS\n"
	 << "#  NVA\n"
         << setw(6) << sys.vas().numVirtualAtoms() << "\n"
         << "# ATNM   IAC       CG   TVA     NVAD VAD[1..NVAD] \n"
         << "#                                INE JNE[1..INE]\n"
         << "#                              INE14 JNE14[1..INE14]\n";
    num = sys.vas().numVirtualAtoms();
    int offatom=0;
    for(unsigned int i=0; i< sys.numMolecules(); i++) 
      offatom += sys.mol(i).numAtoms();
    for(unsigned int i=0; i < num; i++){
      d_os << setw(6) << offatom + i +1
           << setw(6) << sys.vas().iac(i) + 1
           << setw(9) << sys.vas().charge(i)
           << setw(6) << sys.vas().atom(i).type()  
           << setw(9) << sys.vas().atom(i).conf().size();
      for(unsigned int j=0; j< sys.vas().atom(i).conf().size(); j++){
        if (j % 8 == 0 && j != 0)
          d_os << "\n"
                << "                                   ";

        d_os << " " << setw(4) << sys.vas().atom(i).conf().gromosAtom(j) +1 ;
      }
      d_os << "\n";
      gcore::Exclusion e=sys.vas().exclusion(i);
      d_os << "                              " << setw(6) << e.size();
      for(unsigned int j=0; j < e.size(); j++){
        if(j %8 == 0 && j!= 0)
          d_os << "\n"
                << "                                   ";
        d_os << " " << setw(4) << e.atom(j)+1;
      }
      d_os << "\n";
      e=sys.vas().exclusion14(i);
      d_os << "                              " << setw(6) << e.size();
      for(unsigned int j=0; j < e.size(); j++){
        if(j %8 == 0 && j!= 0)
          d_os << "\n"
                << "                                   ";
        d_os << " " << setw(4) << e.atom(j)+1;
      }
      d_os << "\n";


    }
    d_os << "END\n";
  }
  if (args::Arguments::outG96) {
    // BONDTYPE block
    d_os << "BONDTYPE\n"
            << "#  NBTY: number of covalent bond types\n";
    num = gff.numBondTypes();

    d_os << num << "\n"
            << "#  CB: force constant\n"
            << "#  B0: bond length at minimum energy\n"
            << "#         CB          B0\n";

    for (int i = 0; i < num; ++i) {
      if (i > 0 && !(i % 10)) d_os << "# " << i << "\n";
      d_os.precision(5);
      d_os.setf(ios::fixed, ios::floatfield);
      d_os.setf(ios::scientific, ios::floatfield);
      d_os << setw(12) << gff.bondType(i).fc()
           << ' ' << setw(11) << gff.bondType(i).b0() << "\n";
    }
    d_os << "END\n";
  } else {
    // BONDSTRETCHTYPE block

    d_os << "BONDSTRETCHTYPE\n"
            << "#  NBTY: number of covalent bond types\n";
    num = gff.numBondTypes();

    d_os << num << "\n"
            << "#  CB:  quartic force constant\n"
            << "#  CHB: harmonic force constant\n"
            << "#  B0:  bond length at minimum energy\n"
            << "#         CB         CHB         B0\n";

    for (int i = 0; i < num; ++i) {
      if (i > 0 && !(i % 10)) d_os << "# " << i << "\n";
      d_os.precision(5);
      d_os.setf(ios::fixed, ios::floatfield);
      d_os.setf(ios::scientific, ios::floatfield);
      d_os << setw(16) << gff.bondType(i).fc()
           << ' ' << setw(15) << gff.bondType(i).hfc()
           << ' ' << setw(15) << gff.bondType(i).b0() << "\n";
    }
    d_os << "END\n";
  }

  // BONDH and BOND blocks
  num = 0;
  for (int i = 0; i < sys.numMolecules(); ++i) {
    BondIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
      if ((sys.mol(i).topology().atom(bit()[0]).isH() ||
           sys.mol(i).topology().atom(bit()[1]).isH()) &&
         (!sys.mol(i).topology().atom(bit()[0]).isCoarseGrained() &&
          !sys.mol(i).topology().atom(bit()[1]).isCoarseGrained()))
        ++num;
    }
  }
  d_os << "BONDH\n"
          << "#  NBONH: number of bonds involving H atoms in solute\n"
          << num << "\n"
          << "#  IBH, JBH: atom sequence numbers of atoms forming a bond\n"
          << "#  ICBH: bond type code\n"
          << "#   IBH    JBH ICBH\n";

  for (int i = 0, offatom = 1, count = 0; i < sys.numMolecules(); ++i) {
    BondIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
      if ((sys.mol(i).topology().atom(bit()[0]).isH() ||
           sys.mol(i).topology().atom(bit()[1]).isH()) &&
         (!sys.mol(i).topology().atom(bit()[0]).isCoarseGrained() &&
          !sys.mol(i).topology().atom(bit()[1]).isCoarseGrained())) {
        if (count > 0 && !(count % 10))d_os << "# " << count << "\n";
        d_os << setw(7) << bit()[0] + offatom
             << ' ' << setw(6) << bit()[1] + offatom
             << ' ' << setw(4) << bit().type() + 1 << "\n";
        ++count;
      }
    }
    offatom += sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  d_os << "BOND\n"
          << "#  NBON: number of bonds NOT involving H atoms in solute\n";
  num = 0;
  for (int i = 0; i < sys.numMolecules(); ++i) {
    BondIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
      if (!sys.mol(i).topology().atom(bit()[0]).isH() &&
          !sys.mol(i).topology().atom(bit()[1]).isH())
        ++num;
    }
  }
  d_os << num << "\n"
          << "#  IB, JB: atom sequence numbers of atoms forming a bond\n"
          << "#  ICB: bond type code\n"
          << "#    IB     JB  ICB\n";
  
  for (int i = 0, offatom = 1, count = 0; i < sys.numMolecules(); ++i) {
    BondIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
      if (!sys.mol(i).topology().atom(bit()[0]).isH() &&
          !sys.mol(i).topology().atom(bit()[1]).isH()) {
        if (count > 0 && !(count % 10))d_os << "# " << count << "\n";
        d_os << setw(7) << bit()[0] + offatom
             << ' ' << setw(6) << bit()[1] + offatom
             << ' ' << setw(4) << bit().type() + 1 << "\n";
        ++count;
      }
    }
    offatom += sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  // CONSTRAINT block
  num = 0;
  for (int m = 0; m < sys.numMolecules(); ++m) {
    num+=sys.mol(m).topology().constraints().size();
  }
  if (num > 0) {
    d_os << "CONSTRAINT\n"
          << "#  NCON: number of constraints in solute\n"
          << num << "\n"
          << "#  IC, JC: atom sequence numbers of atoms defining the constraint\n"
          << "#  ICC: bond type code\n"
          << "#    IC     JC  ICC\n";
 
    for (int i = 0, offatom = 1, count = 0; i < sys.numMolecules(); ++i) {
      set<Constraint>::const_iterator iter = sys.mol(i).topology().constraints().begin(),
                                    to = sys.mol(i).topology().constraints().end();
      for (; iter != to; ++iter) {
        if (count > 0 && !(count % 10))d_os << "# " << count << "\n";
        d_os << setw(7) << (*iter)[0] + offatom
             << ' ' << setw(6) << (*iter)[1] + offatom
             << ' ' << setw(4) << iter->bondtype() + 1 << "\n";
        ++count;
      }
      offatom += sys.mol(i).numAtoms();
    }
    d_os << "END\n";
  }

  // BONDDP
  //check whether we have CG bonds
  num = 0;
  for (int i = 0; i < sys.numMolecules(); ++i) {
    BondDipoleIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
//      if (sys.mol(i).topology().atom(bit()[0]).isCoarseGrained() &&
//          sys.mol(i).topology().atom(bit()[1]).isCoarseGrained())
        ++num;
    }
  }
  if (num) {
    d_os << "BONDDP\n"
            << "#  NBONDP: number of bonds involving dipole particles in solute\n";
    d_os << num << "\n"
            << "#  IBDP, JBDP: sequence numbers of dipole particles forming a bond\n"
            << "#  ICB: bond type code\n"
            << "#    IBDP   JBDP  ICB\n";

    for (int i = 0, offatom = 1, count = 0; i < sys.numMolecules(); ++i) {
      BondDipoleIterator bit(sys.mol(i).topology());
      for (; bit; ++bit) {
//        if (sys.mol(i).topology().atom(bit()[0]).isCoarseGrained() &&
//            sys.mol(i).topology().atom(bit()[1]).isCoarseGrained()) {
          if (count > 0 && !(count % 10))d_os << "# " << count << "\n";
          d_os << setw(7) << bit()[0] + offatom
               << ' ' << setw(6) << bit()[1] + offatom
               << ' ' << setw(4) << bit().type() + 1 << "\n";
          ++count;
//        }
      }
      offatom += sys.mol(i).numAtoms();
    }
    d_os << "END\n";
  } // BONDDP
  
  if (args::Arguments::outG96) {
    // BONDANGLETYPE block
    num = gff.numAngleTypes();
    d_os << "BONDANGLETYPE\n"
            << "#  NTTY: number of bond angle types\n"
            << num << "\n"
            << "#  CT: force constant\n"
            << "#  T0: bond angle at minimum energy in degrees\n"
            << "#         CT          T0\n";

    for (int i = 0; i < num; ++i) {
      if (i > 0 && !(i % 10))d_os << "# " << i << "\n";
      d_os.precision(5);
      d_os.setf(ios::fixed, ios::floatfield);
      d_os.setf(ios::scientific, ios::floatfield);
      d_os << setw(12) << gff.angleType(i).fc()
           << ' ' << setw(11) << gff.angleType(i).t0() << "\n";
    }
    d_os << "END\n";
  } else {
    // BONDANGLEBENDTYPE block
    num = gff.numAngleTypes();
    d_os << "BONDANGLEBENDTYPE\n"
            << "#  NTTY: number of bond angle types\n"
            << num << "\n"
            << "#  CT:  force constant (based on potential\n"
            << "#       harmonic in the angle cosine)\n"
            << "#  CHT: force constant (based on potential\n"
            << "#       harmonic in the angle)\n"
            << "#  T0:  bond angle at minimum energy in degrees\n"
            << "#         CT         CHT          T0\n";

    for (int i = 0; i < num; ++i) {
      if (i > 0 && !(i % 10))d_os << "# " << i << "\n";
      d_os.precision(5);
      d_os.setf(ios::fixed, ios::floatfield);
      d_os.setf(ios::scientific, ios::floatfield);
      d_os << setw(16) << gff.angleType(i).fc()
           << ' ' << setw(15) << gff.angleType(i).afc()
           << ' ' << setw(15) << gff.angleType(i).t0() << "\n";
    }
    d_os << "END\n";
  }

  // BONDANGLEH & BONDANGLE block
  num = 0;
  for (int i = 0; i < sys.numMolecules(); ++i) {
    AngleIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
      if (sys.mol(i).topology().atom(bit()[0]).isH() ||
              sys.mol(i).topology().atom(bit()[1]).isH() ||
              sys.mol(i).topology().atom(bit()[2]).isH())
        ++num;
    }
  } 
  d_os << "BONDANGLEH\n"
          << "#  NTHEH: number of bond angles involving H atoms in solute\n"
          << num << "\n"
          << "#  ITH, JTH, KTH: atom sequence numbers\n"
          << "#    of atoms forming a bond angle in solute\n"
          << "#  ICTH: bond angle type code\n"
          << "#   ITH    JTH    KTH ICTH\n";

  for (int i = 0, offatom = 1, count = 0; i < sys.numMolecules(); ++i) {
    AngleIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
      if (sys.mol(i).topology().atom(bit()[0]).isH() ||
              sys.mol(i).topology().atom(bit()[1]).isH() ||
              sys.mol(i).topology().atom(bit()[2]).isH()) {

        if (count > 0 && !(count % 10))d_os << "# " << count << "\n";
        d_os << setw(7) << bit()[0] + offatom
             << ' ' << setw(6) << bit()[1] + offatom
             << ' ' << setw(6) << bit()[2] + offatom
             << ' ' << setw(4) << bit().type() + 1 << "\n";
        ++count;
      }
    }
    offatom += sys.mol(i).numAtoms();
  }
  d_os << "END\n";
  
  num = 0;
  for (int i = 0; i < sys.numMolecules(); ++i) {
    AngleIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
      if (!sys.mol(i).topology().atom(bit()[0]).isH() &&
              !sys.mol(i).topology().atom(bit()[1]).isH() &&
              !sys.mol(i).topology().atom(bit()[2]).isH())
        ++num;
    }
  }
  d_os << "BONDANGLE\n"
          << "#  NTHE: number of bond angles NOT\n"
          << "#        involving H atoms in solute\n"
          << num << "\n"
          << "#  IT, JT, KT: atom sequence numbers of atoms\n"
          << "#     forming a bond angle\n"
          << "#  ICT: bond angle type code\n"
          << "#    IT     JT     KT  ICT\n";

  for (int i = 0, offatom = 1, count = 0; i < sys.numMolecules(); ++i) {
    AngleIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
      if (!sys.mol(i).topology().atom(bit()[0]).isH() &&
              !sys.mol(i).topology().atom(bit()[1]).isH() &&
              !sys.mol(i).topology().atom(bit()[2]).isH()) {
        if (count > 0 && !(count % 10))d_os << "# " << count << "\n";
        d_os << setw(7) << bit()[0] + offatom
             << ' ' << setw(6) << bit()[1] + offatom
             << ' ' << setw(6) << bit()[2] + offatom
             << ' ' << setw(4) << bit().type() + 1 << "\n";
        ++count;
      }
    }
    offatom += sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  // IMPDIHEDRALTYPE block
  num = gff.numImproperTypes();
  d_os << "IMPDIHEDRALTYPE\n"
          << "#  NQTY: number of improper dihedrals\n"
          << num << "\n"
          << "#  CQ: force constant of improper dihedral per degrees square\n"
          << "#  Q0: improper dihedral angle at minimum energy in degrees\n"
          << "#            CQ             Q0\n";

  for (int i = 0; i < num; ++i) {
    if (i > 0 && !(i % 10))d_os << "# " << i << "\n";
    d_os.precision(5);
    d_os.setf(ios::fixed, ios::floatfield);
    d_os.setf(ios::scientific, ios::floatfield);
    d_os << setw(15) << gff.improperType(i).fc()
         << ' ' << setw(14) << gff.improperType(i).q0() << "\n";
  }
  d_os << "END\n";

  // IMPDIHEDRALH & IMPDIHEDRAL block
  num = 0;
  for (int i = 0; i < sys.numMolecules(); ++i) {
    ImproperIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
      if (sys.mol(i).topology().atom(bit()[0]).isH() ||
              sys.mol(i).topology().atom(bit()[1]).isH() ||
              sys.mol(i).topology().atom(bit()[2]).isH() ||
              sys.mol(i).topology().atom(bit()[3]).isH())
        ++num;
    }
  }
  d_os << "IMPDIHEDRALH\n"
          << "#  NQHIH: number of improper dihedrals\n"
          << "#         involving H atoms in the solute\n"
          << num << "\n"
          << "#  IQH,JQH,KQH,LQH: atom sequence numbers\n"
          << "#     of atoms forming an improper dihedral\n"
          << "#  ICQH: improper dihedral type code\n"
          << "#   IQH    JQH    KQH    LQH ICQH\n";

  for (int i = 0, offatom = 1, count = 0; i < sys.numMolecules(); ++i) {
    ImproperIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
      if (sys.mol(i).topology().atom(bit()[0]).isH() ||
              sys.mol(i).topology().atom(bit()[1]).isH() ||
              sys.mol(i).topology().atom(bit()[2]).isH() ||
              sys.mol(i).topology().atom(bit()[3]).isH()) {
        if (count > 0 && !(count % 10))d_os << "# " << count << "\n";
        d_os << setw(7) << bit()[0] + offatom
             << ' ' << setw(6) << bit()[1] + offatom
             << ' ' << setw(6) << bit()[2] + offatom
             << ' ' << setw(6) << bit()[3] + offatom
             << ' ' << setw(4) << bit().type() + 1 << "\n";
        ++count;
      }
    }
    offatom += sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  num = 0;
  for (int i = 0; i < sys.numMolecules(); ++i) {
    ImproperIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
      if (!sys.mol(i).topology().atom(bit()[0]).isH() &&
              !sys.mol(i).topology().atom(bit()[1]).isH() &&
              !sys.mol(i).topology().atom(bit()[2]).isH() &&
              !sys.mol(i).topology().atom(bit()[3]).isH())
        ++num;
    }
  }
  d_os << "IMPDIHEDRAL\n"
          << "#  NQHI: number of improper dihedrals NOT\n"
          << "#    involving H atoms in solute\n"
          << num << "\n"
          << "#  IQ,JQ,KQ,LQ: atom sequence numbers of atoms\n"
          << "#    forming an improper dihedral\n"
          << "#  ICQ: improper dihedral type code\n"
          << "#    IQ     JQ     KQ     LQ  ICQ\n";

  for (int i = 0, offatom = 1, count = 0; i < sys.numMolecules(); ++i) {
    ImproperIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
      if (!sys.mol(i).topology().atom(bit()[0]).isH() &&
              !sys.mol(i).topology().atom(bit()[1]).isH() &&
              !sys.mol(i).topology().atom(bit()[2]).isH() &&
              !sys.mol(i).topology().atom(bit()[3]).isH()) {
        if (count > 0 && !(count % 10))d_os << "# " << count << "\n";
        d_os << setw(7) << bit()[0] + offatom
             << ' ' << setw(6) << bit()[1] + offatom
             << ' ' << setw(6) << bit()[2] + offatom
             << ' ' << setw(6) << bit()[3] + offatom
             << ' ' << setw(4) << bit().type() + 1 << "\n";
        ++count;
      }
    }
    offatom += sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  if (args::Arguments::outG96) {
    // DIHEDRALTYPE block
    num = gff.numDihedralTypes();

    d_os << "DIHEDRALTYPE\n"
            << "#  NPTY: number of dihedral types\n"
            << num << "\n"
            << "#  CP: force constant\n"
            << "#  PD: cosine of the phase shift\n"
            << "#  NP: multiplicity\n"
            << "#       CP        PD  NP\n";

    for (int i = 0; i < num; ++i) {
      if (i > 0 && !(i % 10))d_os << "# " << i << "\n";
      d_os.precision(5);
      d_os.setf(ios::fixed, ios::floatfield);
      d_os << setw(10) << gff.dihedralType(i).fc()
           << ' ' << setw(9) << gff.dihedralType(i).pd()
           << ' ' << setw(3) << gff.dihedralType(i).np() << "\n";
    }
    d_os << "END\n";
  } else {
    // TORSDIHEDRALTYPE block
    num = gff.numDihedralTypes();

    d_os << "TORSDIHEDRALTYPE\n"
            << "#  NPTY: number of dihedral types\n"
            << num << "\n"
            << "#  CP: force constant\n"
            << "#  PD: phase-shift angle\n"
            << "#  NP: multiplicity\n"
            << "#       CP         PD  NP\n";

    for (int i = 0; i < num; ++i) {
      if (i > 0 && !(i % 10)) d_os << "# " << i << "\n";
      d_os.precision(5);
      d_os.setf(ios::fixed, ios::floatfield);
      d_os << setw(10) << gff.dihedralType(i).fc()
           << " " << setw(10) << gff.dihedralType(i).pdl()
           << ' ' << setw(3) << gff.dihedralType(i).np() << "\n";
    }
    d_os << "END\n";
  }

  // DIHEDRALH & DIHEDRAL block
  num = 0;
  for (int i = 0; i < sys.numMolecules(); ++i) {
    DihedralIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
      if (sys.mol(i).topology().atom(bit()[0]).isH() ||
              sys.mol(i).topology().atom(bit()[1]).isH() ||
              sys.mol(i).topology().atom(bit()[2]).isH() ||
              sys.mol(i).topology().atom(bit()[3]).isH())
        ++num;
    }
  }
  d_os << "DIHEDRALH\n"
          << "#  NPHIH: number of dihedrals involving H atoms in solute\n"
          << num << "\n"
          << "#  IPH, JPH, KPH, LPH: atom sequence numbers\n"
          << "#    of atoms forming a dihedral\n"
          << "#  ICPH: dihedral type code\n"
          << "#   IPH    JPH    KPH    LPH ICPH\n";

  for (int i = 0, offatom = 1, count = 0; i < sys.numMolecules(); ++i) {
    DihedralIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
      if (sys.mol(i).topology().atom(bit()[0]).isH() ||
              sys.mol(i).topology().atom(bit()[1]).isH() ||
              sys.mol(i).topology().atom(bit()[2]).isH() ||
              sys.mol(i).topology().atom(bit()[3]).isH()) {
        if (count > 0 && !(count % 10))d_os << "# " << count << "\n";
        d_os << setw(7) << bit()[0] + offatom
             << ' ' << setw(6) << bit()[1] + offatom
             << ' ' << setw(6) << bit()[2] + offatom
             << ' ' << setw(6) << bit()[3] + offatom
             << ' ' << setw(4) << bit().type() + 1 << "\n";
        ++count;
      }
    }
    offatom += sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  num = 0;
  for (int i = 0; i < sys.numMolecules(); ++i) {
    DihedralIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
      if (!sys.mol(i).topology().atom(bit()[0]).isH() &&
              !sys.mol(i).topology().atom(bit()[1]).isH() &&
              !sys.mol(i).topology().atom(bit()[2]).isH() &&
              !sys.mol(i).topology().atom(bit()[3]).isH())
        ++num;
    }
  }
  d_os << "DIHEDRAL\n"
          << "#  NPHI: number of dihedrals NOT involving H atoms in solute\n"
          << num << "\n"
          << "#  IP, JP, KP, LP: atom sequence numbers\n"
          << "#     of atoms forming a dihedral\n"
          << "#  ICP: dihedral type code\n"
          << "#    IP     JP     KP     LP  ICP\n";

  for (int i = 0, offatom = 1, count = 0; i < sys.numMolecules(); ++i) {
    DihedralIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
      if (!sys.mol(i).topology().atom(bit()[0]).isH() &&
              !sys.mol(i).topology().atom(bit()[1]).isH() &&
              !sys.mol(i).topology().atom(bit()[2]).isH() &&
              !sys.mol(i).topology().atom(bit()[3]).isH()) {
        if (count > 0 && !(count % 10))d_os << "# " << count << "\n";
        d_os << setw(7) << bit()[0] + offatom
             << ' ' << setw(6) << bit()[1] + offatom
             << ' ' << setw(6) << bit()[2] + offatom
             << ' ' << setw(6) << bit()[3] + offatom
             << ' ' << setw(4) << bit().type() + 1 << "\n";
        ++count;
      }
    }
    offatom += sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  // CROSSDIHEDRALH & CROSSDIHEDRAL block
  if (!args::Arguments::outG96) {
    num = 0;
    for (int i = 0; i < sys.numMolecules(); ++i) {
      CrossDihedralIterator bit(sys.mol(i).topology());
      for (; bit; ++bit) {
        if (sys.mol(i).topology().atom(bit()[0]).isH() ||
                sys.mol(i).topology().atom(bit()[1]).isH() ||
                sys.mol(i).topology().atom(bit()[2]).isH() ||
                sys.mol(i).topology().atom(bit()[3]).isH() ||
                sys.mol(i).topology().atom(bit()[4]).isH() ||
                sys.mol(i).topology().atom(bit()[5]).isH() ||
                sys.mol(i).topology().atom(bit()[6]).isH() ||
                sys.mol(i).topology().atom(bit()[7]).isH())
          ++num;
      }
    }
    d_os << "CROSSDIHEDRALH\n"
            << "#  NPHIH: number of cross dihedrals involving H atoms in solute\n"
            << num << "\n"
            << "#  APH, BPH, CPH, DPH, EPH, FPH, GPH, HPH: atom sequence numbers\n"
            << "#    of atoms forming a dihedral\n"
            << "#  ICCH: dihedral type code\n"
            << "#   APH    BPH    CPH    DPH    EPH    FPH    GPH    HPH ICCH\n";

    for (int i = 0, offatom = 1, count = 0; i < sys.numMolecules(); ++i) {
      CrossDihedralIterator bit(sys.mol(i).topology());
      for (; bit; ++bit) {
        if (sys.mol(i).topology().atom(bit()[0]).isH() ||
                sys.mol(i).topology().atom(bit()[1]).isH() ||
                sys.mol(i).topology().atom(bit()[2]).isH() ||
                sys.mol(i).topology().atom(bit()[3]).isH() ||
                sys.mol(i).topology().atom(bit()[4]).isH() ||
                sys.mol(i).topology().atom(bit()[5]).isH() ||
                sys.mol(i).topology().atom(bit()[6]).isH() ||
                sys.mol(i).topology().atom(bit()[7]).isH()) {
          if (count > 0 && !(count % 10))d_os << "# " << count << "\n";
          d_os << setw(7) << bit()[0] + offatom
                  << setw(7) << bit()[1] + offatom
                  << setw(7) << bit()[2] + offatom
                  << setw(7) << bit()[3] + offatom
                  << setw(7) << bit()[4] + offatom
                  << setw(7) << bit()[5] + offatom
                  << setw(7) << bit()[6] + offatom
                  << setw(7) << bit()[7] + offatom
                  << setw(5) << bit().type() + 1 << "\n";
          ++count;
        }
      }
      offatom += sys.mol(i).numAtoms();
    }
    d_os << "END\n";

    num = 0;
    for (int i = 0; i < sys.numMolecules(); ++i) {
      CrossDihedralIterator bit(sys.mol(i).topology());
      for (; bit; ++bit) {
        if (!sys.mol(i).topology().atom(bit()[0]).isH() &&
                !sys.mol(i).topology().atom(bit()[1]).isH() &&
                !sys.mol(i).topology().atom(bit()[2]).isH() &&
                !sys.mol(i).topology().atom(bit()[3]).isH() &&
                !sys.mol(i).topology().atom(bit()[4]).isH() &&
                !sys.mol(i).topology().atom(bit()[5]).isH() &&
                !sys.mol(i).topology().atom(bit()[6]).isH() &&
                !sys.mol(i).topology().atom(bit()[7]).isH())
          ++num;
      }
    }
    d_os << "CROSSDIHEDRAL\n"
            << "#  NPPC: number of cross dihedrals NOT involving H atoms in solute\n"
            << num << "\n"
            << "#  AP, BP, CP, DP, EP, FP, GP, HP: atom sequence numbers\n"
            << "#     of atoms forming a dihedral\n"
            << "#  ICC: dihedral type code\n"
            << "#    AP     BP     CP     DP     EP     FP     GP     HP  ICC\n";

    for (int i = 0, offatom = 1, count = 0; i < sys.numMolecules(); ++i) {
      CrossDihedralIterator bit(sys.mol(i).topology());
      for (; bit; ++bit) {
        if (!sys.mol(i).topology().atom(bit()[0]).isH() &&
                !sys.mol(i).topology().atom(bit()[1]).isH() &&
                !sys.mol(i).topology().atom(bit()[2]).isH() &&
                !sys.mol(i).topology().atom(bit()[3]).isH() &&
                !sys.mol(i).topology().atom(bit()[4]).isH() &&
                !sys.mol(i).topology().atom(bit()[5]).isH() &&
                !sys.mol(i).topology().atom(bit()[6]).isH() &&
                !sys.mol(i).topology().atom(bit()[7]).isH()) {
          if (count > 0 && !(count % 10))d_os << "# " << count << "\n";
          d_os << setw(7) << bit()[0] + offatom
                  << setw(7) << bit()[1] + offatom
                  << setw(7) << bit()[2] + offatom
                  << setw(7) << bit()[3] + offatom
                  << setw(7) << bit()[4] + offatom
                  << setw(7) << bit()[5] + offatom
                  << setw(7) << bit()[6] + offatom
                  << setw(7) << bit()[7] + offatom
                  << setw(5) << bit().type() + 1 << "\n";
          ++count;
        }
      }
      offatom += sys.mol(i).numAtoms();
    }
    d_os << "END\n";
  } else {
    num = 0;
    for (int i = 0; i < sys.numMolecules(); ++i) {
      CrossDihedralIterator bit(sys.mol(i).topology());
      for (; bit; ++bit) {
        if (sys.mol(i).topology().atom(bit()[0]).isH() ||
                sys.mol(i).topology().atom(bit()[1]).isH() ||
                sys.mol(i).topology().atom(bit()[2]).isH() ||
                sys.mol(i).topology().atom(bit()[3]).isH() ||
                sys.mol(i).topology().atom(bit()[4]).isH() ||
                sys.mol(i).topology().atom(bit()[5]).isH() ||
                sys.mol(i).topology().atom(bit()[6]).isH() ||
                sys.mol(i).topology().atom(bit()[7]).isH())
          ++num;
      }
    }
    for (int i = 0; i < sys.numMolecules(); ++i) {
      CrossDihedralIterator bit(sys.mol(i).topology());
      for (; bit; ++bit) {
        if (!sys.mol(i).topology().atom(bit()[0]).isH() &&
                !sys.mol(i).topology().atom(bit()[1]).isH() &&
                !sys.mol(i).topology().atom(bit()[2]).isH() &&
                !sys.mol(i).topology().atom(bit()[3]).isH() &&
                !sys.mol(i).topology().atom(bit()[4]).isH() &&
                !sys.mol(i).topology().atom(bit()[5]).isH() &&
                !sys.mol(i).topology().atom(bit()[6]).isH() &&
                !sys.mol(i).topology().atom(bit()[7]).isH())
          ++num;
      }
    }
    if (num > 0) {
      stringstream msg;
      msg << "cross dihedral(s) present but not supported in the old format (outG96)";
      throw gromos::Exception("OutTopology", msg.str());
    }
  }

  // LJPARAMETERS block
  num = gff.numLJTypes();
  d_os << "LJPARAMETERS\n"
          << "#  NRATT2: number of LJ interaction types = NRATT*(NRATT+1)/2\n"
          << num << "\n"
          << "#  IAC,JAC: integer (van der Waals) atom type code\n"
          << "#  C12: r**(-12) term in nonbonded interactions\n"
          << "#   C6: r**(-6) term in nonbonded interactions\n"
          << "# CS12: r**(-12) term in 1-4 nonbonded interactions\n"
          << "#  CS6: r**(-6) term in 1-4 nonbonded interactions\n"
          << "# IAC  JAC           C12            C6          CS12           CS6\n";

  for (int i = 0; i < gff.numAtomTypeNames(); ++i) {
    for (int j = 0; j <= i; ++j) {
      d_os.precision(6);
      d_os.setf(ios::fixed, ios::floatfield);
      d_os.setf(ios::scientific, ios::floatfield);
      LJType lj(gff.ljType(AtomPair(i, j)));
      d_os << setw(5) << j + 1
              << setw(5) << i + 1
              << setw(14) << lj.c12()
              << setw(14) << lj.c6()
              << setw(14) << lj.cs12()
              << setw(14) << lj.cs6() << "\n";
    }
    d_os << "#\n";
  }
  d_os << "END\n";

  // CGPARAMETERS block
  num = gff.numCGTypes();
  // only write the block if there are parameters and not in G96 mode.
  if (num && !(args::Arguments::outG96)) {
    d_os << "CGPARAMETERS\n"
            << "#  NRATT2: number of coarse grain LJ interaction types = NRATT*(NRATT+1)/2\n"
            << num << "\n"
            << "#  IAC,JAC: integer (van der Waals) atom type code\n"
            << "#  C12: r**(-12) term in nonbonded interactions\n"
            << "#   C6: r**(-6) term in nonbonded interactions\n"
            << "# IAC  JAC           C12            C6\n";

    for (int i = 0; i < gff.numAtomTypeNames(); ++i) {
      for (int j = 0; j <= i; ++j) {
        d_os.precision(6);
        d_os.setf(ios::fixed, ios::floatfield);
        d_os.setf(ios::scientific, ios::floatfield);
        CGType cg(gff.cgType(AtomPair(i, j)));
        d_os << setw(5) << j + 1
                << setw(5) << i + 1
                << setw(14) << cg.c12()
                << setw(14) << cg.c6() << "\n";
      }
      d_os << "#\n";
    }
    d_os << "END\n";
  }


  if (!(args::Arguments::outG96)) {
    // SOLUTEMOLECULES block
    // Still implemented for default only
    d_os << "SOLUTEMOLECULES\n"
            << "# NSPM: number of separate molecules in solute block\n"
            << "# NSP[1...NSPM]: atom sequence number of last atom\n"
            << "#                of the successive submolecules\n"
            << "#      NSPM  NSP[1...NSPM]\n";

    d_os << setw(10) << sys.numMolecules() << "\n";
    int nspmin = 0;
    for (int i = 0; i < sys.numMolecules(); ++i) {
      d_os << ' ' << setw(5) << sys.mol(i).numAtoms() + nspmin;
      nspmin += sys.mol(i).numAtoms();
      if ((i + 1) % 10 == 0) d_os << "\n";
    }

    if (sys.numMolecules() % 10 != 0) d_os << "\n";
    d_os << "END\n";

    // TEMPERATUREGROUPS block
    // Still implemented for default only
    d_os << "TEMPERATUREGROUPS\n"
            << "# NSTM: number of temperature atom groups\n"
            << "# NST[1...NSTM]: atom sequence number of last atom\n"
            << "#                of the successive temperature atom groups\n"
            << "#      NSTM  NST[1...NSTM]\n";

    d_os << setw(10) << sys.numTemperatureGroups() << "\n";
    for (int i = 0; i < sys.numTemperatureGroups(); ++i) {
      d_os << ' ' << setw(5) << sys.temperatureGroup(i);
      if ((i + 1) % 10 == 0) d_os << "\n";
    }

    if (sys.numTemperatureGroups() % 10 != 0) d_os << "\n";
    d_os << "END\n";

    // PRESSUREGROUPS block
    // Still implemented for default only
    d_os << "PRESSUREGROUPS\n"
            << "# NSVM: number of pressure atom groups\n"
            << "# NSV[1...NSVM]: atom sequence number of last atom\n"
            << "#                of the successive pressure atom groups\n"
            << "#      NSVM  NSV[1...NSVM]\n";

    d_os << setw(10) << sys.numPressureGroups() << "\n";
    for (int i = 0; i < sys.numPressureGroups(); ++i) {
      d_os << ' ' << setw(5) << sys.pressureGroup(i);
      if ((i + 1) % 10 == 0) d_os << "\n";
    }

    if (sys.numPressureGroups() % 10 != 0) d_os << "\n";
    d_os << "END\n";
  }

  // LJEXCEPTIONS
  if (!args::Arguments::outG96) {
    d_os << "LJEXCEPTIONS\n"
            << "# This block defines special LJ-interactions based on atom numbers\n"
            << "# This overrules the normal LJ-parameters (including 1-4 interactions)\n"
            << "# NEX: number of exceptions\n";
    num = 0;
    for (int i = 0; i < sys.numMolecules(); ++i) {
      LJExceptionIterator ljit(sys.mol(i).topology());
      for (; ljit; ++ljit) {
        ++num;
      }
    }
    d_os << num << "\n"
            << "# AT1  AT2           C12            C6\n";
    for (int i = 0, offatom = 1; i < sys.numMolecules(); ++i) {
      LJExceptionIterator ljit(sys.mol(i).topology());
      for (int count = 0; ljit; ++ljit) {
        if (count > 0 && !(count % 10))d_os << "# " << count << "\n";
        // get the C12 and C6 parameters for these LJ exception (from the types)
        int found = 0;
        int numLJTypes = gff.numLJExceptionTypes();
        double c12, c6;
        for (int j = 0; j < numLJTypes; ++j) {
          if (ljit().type() == gff.ljExceptionType(j).code()) {
            found = 1;
            c12 = gff.ljExceptionType(j).c12();
            c6 = gff.ljExceptionType(j).c6();
            break;
          }
        }
        if (found == 0) {
          stringstream msg;
          msg << "no C12 and C6 parameters for LJ exception of type " << ljit().type() + 1 << " found";
          throw gromos::Exception("OutTopology", msg.str());
        }
        d_os << setw(5) << ljit()[0] + offatom
             << ' ' << setw(4) << ljit()[1] + offatom
             << ' ' << setw(13) << c12 
             << ' ' << setw(13) << c6 << "\n";
        ++count;
      }
      offatom += sys.mol(i).numAtoms();
    }
    d_os << "END\n";
  } else {
    num = 0;
    for (int i = 0; i < sys.numMolecules(); ++i) {
      LJExceptionIterator ljit(sys.mol(i).topology());
      for (; ljit; ++ljit) {
        ++num;
      }
    }
    if (num > 0) {
      stringstream msg;
      msg << "there were LJ exceptions which are not supported in the desired format (outG96)";
      throw gromos::Exception("OutTopology", msg.str());
    }
  }

  //SOLVENTATOM block
  d_os << "SOLVENTATOM\n"
          << "#  NRAM: number of atoms per solvent molecule\n"
          << sys.sol(0).topology().numAtoms() << "\n"
          << "#     I: solvent atom sequence number\n"
          << "#  IACS: integer (van der Waals) atom type code\n"
          << "#  ANMS: atom name of solvent atom\n"
          << "#  MASS: mass of solvent atom\n"
          << "#   CGS: charge of solvent atom\n"
          << "#  I  ANMS IACS      MASS        CGS\n";

  for (int j = 0; j < sys.sol(0).topology().numAtoms(); ++j) {
    d_os.precision(5);
    d_os.setf(ios::fixed, ios::floatfield);
    d_os << setw(4) << j + 1 << " "
         << ' ' << setw(4) << sys.sol(0).topology().atom(j).name()
         << ' ' << setw(3) << sys.sol(0).topology().atom(j).iac() + 1
         << ' ' << setw(10) << sys.sol(0).topology().atom(j).mass()
         << ' ' << setw(10) << sys.sol(0).topology().atom(j).charge();
    d_os << "\n";
  }
  d_os << "END\n";

  // SOLVENTPOLARISATION block
  // check whether we have polarisable solvent atoms;
  num = 0;
  for (int i = 0; i < sys.sol(0).topology().numAtoms(); ++i) {
    if (sys.sol(0).topology().atom(i).isPolarisable()) ++num;
  }
  if (num) {
    d_os << "SOLVENTPOLARISATION\n"
            << "# NVPOL: number of polarisable solvent atoms\n"
            << setw(5) << num << "\n"
            << "# IPOLV: atom sequence number of the polarisable solvent atom\n"
            << "#  ALPV: polarisability of the solvent atom\n"
            << "# QPOLV: size of charge-on-spring connected to polarisable solvent atoms\n"
            << "# ENOTV: damping level for polarisation\n"
            << "#   EPV: damping power for polarisation\n"
            << "#  GAMV: gamma of for the off-site pol. centre construction\n"
            << "#    IV: first atom for the off-site pol. centre construction\n"
            << "#    JV: second atom for the off-site pol. centre construction\n#\n"
            << "# IPOLV" << setw(14) << "ALPV" << setw(15) << "QPOLV" << setw(15)
            << "ENOTV" << setw(15) << "EPV" << setw(15) << "GAMV"
            << setw(15) << "IV" << setw(15) << "JV" << endl;
    for (int j = 0; j < sys.sol(0).topology().numAtoms(); ++j) {
      if (!sys.sol(0).topology().atom(j).isPolarisable()) continue;

      d_os.precision(7);
      d_os.setf(ios::fixed, ios::floatfield);
      d_os << setw(5) << j + 1 << ' '
              << setw(15) << sys.sol(0).topology().atom(j).polarisability()
              << setw(15) << sys.sol(0).topology().atom(j).cosCharge()
              << setw(15) << sys.sol(0).topology().atom(j).dampingLevel()
              << setw(15) << sys.sol(0).topology().atom(j).dampingPower()
              << setw(15) << sys.sol(0).topology().atom(j).poloffsiteGamma()
              << setw(15) << sys.sol(0).topology().atom(j).poloffsiteI() + 1
              << setw(15) << sys.sol(0).topology().atom(j).poloffsiteJ() + 1
              << "\n";
    }
    d_os << "END\n";
  } // SOLVENTPOLARISATION

  //SOLVENTCONSTR bock
  num = 0;
  ConstraintIterator dit(sys.sol(0).topology());
  for (; dit; ++dit) ++num;
  d_os << "SOLVENTCONSTR\n"
          << "#  NCONS: number of constraints\n"
          << num << "\n"
          << "#  ICONS, JCONS: atom sequence numbers forming constraint\n"
          << "#   CONS constraint length\n"
          << "#ICONS JCONS         CONS\n";

  ConstraintIterator cit(sys.sol(0).topology());

  for (; cit; ++cit) {
    d_os.precision(7);
    d_os.setf(ios::fixed, ios::floatfield);

    d_os << setw(5) << cit()[0] + 1
         << ' ' << setw(4) << cit()[1] + 1
         << ' ' << setw(14) << cit().dist() << "\n";
  }
  d_os << "END\n";
  
  d_os << "# end of topology file" << endl;
}

//OutTopology &OutTopology::operator<<(const gcore::Simulation &sim){

void OutTopology::write96(const gcore::System &sys, const gcore::GromosForceField & gff) {

  // Title block
  d_os << "TITLE\n" << d_title << "\nEND\n";

  // TOPPHYSCON block
  d_os.precision(10);

  d_os << "TOPPHYSCON\n"
          << "# FPEPSI: 1.0/(4.0*PI*EPS0) (EPS0 is the permittivity of vacuum)\n"
          << gff.fpepsi()
          << "\n# HBAR: Planck's constant HBAR = H/(2* PI)\n"
          << gff.hbar()
          << "\n# SPDL: Speed of light (nm/ps)\n"
          << gff.spdl()
          << "\nEND\n";

  // TOPVERSION block
  d_os << "TOPVERSION\n1.7\nEND\n";

  // ATOMTYPENAME block
  d_os << "ATOMTYPENAME\n"
          << "# NRATT: number of van der Waals atom types\n";
  int num = gff.numAtomTypeNames();
  d_os << num << "\n";
  d_os << "# TYPE: atom type names\n";
  for (int i = 0; i < num; ++i) {
    if (i > 0 && !(i % 10))d_os << "# " << i << "\n";
    d_os << gff.atomTypeName(i) << "\n";
  }
  d_os << "END\n";

  // RESNAME block
  d_os << "RESNAME\n"
          << "# NRAA2: number of residues in a solute molecule\n";
  num = 0;
  for (int i = 0; i < sys.numMolecules(); ++i)
    num += sys.mol(i).topology().numRes();

  d_os << num << "\n"
          << "# AANM: residue names\n";
  for (int i = 0, count = 0; i < sys.numMolecules(); ++i)
    for (int j = 0; j < sys.mol(i).topology().numRes(); ++j, ++count) {
      if (count > 0 && !(count % 10))d_os << "# " << count << "\n";
      d_os << sys.mol(i).topology().resName(j) << "\n";
    }
  d_os << "END\n";

  // SOLUTEATOM block
  d_os << "SOLUTEATOM\n"
          << "#   NRP: number of solute atoms\n";
  num = 0;
  for (int i = 0; i < sys.numMolecules(); ++i)
    num += sys.mol(i).numAtoms();

  d_os << setw(5) << num << "\n";
  d_os << "#  ATNM: atom number\n"
          << "#  MRES: residue number\n"
          << "#  PANM: atom name of solute atom\n"
          << "#   IAC: integer (van der Waals) atom type code\n"
          << "#  MASS: mass of solute atom\n"
          << "#    CG: charge of solute atom\n"
          << "#   CGC: charge group code (0 or 1)\n"
          << "#   INE: number of excluded atoms\n"
          << "# INE14: number of 1-4 interactions\n"
          << "# ATNM MRES PANM IAC     MASS       CG  CGC INE\n"
          << "#                                           INE14\n";

  for (int i = 0, offatom = 1, offres = 1; i < sys.numMolecules(); ++i) {
    for (int j = 0; j < sys.mol(i).numAtoms(); ++j) {
      d_os.precision(5);
      d_os.setf(ios::fixed, ios::floatfield);
      d_os << setw(6) << offatom + j << ' '
              << setw(4) << sys.mol(i).topology().resNum(j) + offres << ' '
              << setw(4) << sys.mol(i).topology().atom(j).name()
              << setw(4) << sys.mol(i).topology().atom(j).iac() + 1
              << setw(9) << sys.mol(i).topology().atom(j).mass()
              << setw(9) << sys.mol(i).topology().atom(j).charge()
              << setw(3) << sys.mol(i).topology().atom(j).chargeGroup();
      // Exclusions
      d_os << setw(3) << sys.mol(i).topology().atom(j).exclusion().size();
      for (int k = 0; k < sys.mol(i).topology().atom(j).exclusion().size(); ++k) {
        if (k % 6 == 0 && k != 0)
          d_os << "\n"
                << "                                            ";
        d_os << setw(6) << sys.mol(i).topology().atom(j).exclusion().atom(k) + offatom;
      }
      d_os << "\n"
              << "                                           "
              << sys.mol(i).topology().atom(j).exclusion14().size();
      for (int k = 0; k < sys.mol(i).topology().atom(j).exclusion14().size(); ++k) {
        if (k % 6 == 0 && k != 0)
          d_os << "\n"
                << "                                            ";
        d_os << setw(6) << sys.mol(i).topology().atom(j).exclusion14().atom(k) + offatom;
      }

      d_os << "\n";
    }
    offres += sys.mol(i).topology().numRes();
    offatom += sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  // BONDTYPE block

  d_os << "BONDTYPE\n"
          << "#  NBTY: number of covalent bond types\n";
  num = gff.numBondTypes();

  d_os << num << "\n"
          << "#  CB: force constant\n"
          << "#  B0: bond length at minimum energy\n"
          << "#         CB          B0\n";

  for (int i = 0; i < num; ++i) {
    if (i > 0 && !(i % 10)) d_os << "# " << i << "\n";
    d_os.precision(5);
    d_os.setf(ios::fixed, ios::floatfield);
    d_os.setf(ios::scientific, ios::floatfield);
    d_os << setw(12) << gff.bondType(i).fc()
            << setw(12) << gff.bondType(i).b0() << "\n";
  }
  d_os << "END\n";

  // BONDH block
  num = 0;
  for (int i = 0; i < sys.numMolecules(); ++i) {
    BondIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
      if (sys.mol(i).topology().atom(bit()[0]).isH() ||
              sys.mol(i).topology().atom(bit()[1]).isH())
        ++num;
    }
  }
  d_os << "BONDH\n"
          << "#  NBONH: number of bonds involving H atoms in solute\n"
          << num << "\n"
          << "#  IBH, JBH: atom sequence numbers of atoms forming a bond\n"
          << "#  ICBH: bond type code\n"
          << "#   IBH    JBH ICBH\n";

  for (int i = 0, offatom = 1; i < sys.numMolecules(); ++i) {
    BondIterator bit(sys.mol(i).topology());
    for (int count = 0; bit; ++bit) {
      if (sys.mol(i).topology().atom(bit()[0]).isH() ||
              sys.mol(i).topology().atom(bit()[1]).isH()) {
        if (count > 0 && !(count % 10))d_os << "# " << count << "\n";
        d_os << setw(7) << bit()[0] + offatom
                << setw(7) << bit()[1] + offatom
                << setw(5) << bit().type() + 1 << "\n";
        ++count;
      }
    }
    offatom += sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  d_os << "BOND\n"
          << "#  NBON: number of bonds NOT involving H atoms in solute\n";
  num = 0;
  for (int i = 0; i < sys.numMolecules(); ++i) {
    BondIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
      if (!sys.mol(i).topology().atom(bit()[0]).isH() &&
          !sys.mol(i).topology().atom(bit()[1]).isH())
        ++num;
    }
  }
  d_os << num << "\n"
          << "#  IB, JB: atom sequence numbers of atoms forming a bond\n"
          << "#  ICB: bond type code\n"
          << "#    IB     JB  ICB\n";

  for (int i = 0, offatom = 1; i < sys.numMolecules(); ++i) {
    BondIterator bit(sys.mol(i).topology());
    for (int count = 0; bit; ++bit) {
      if (!sys.mol(i).topology().atom(bit()[0]).isH() &&
          !sys.mol(i).topology().atom(bit()[1]).isH()) {
        if (count > 0 && !(count % 10))d_os << "# " << count << "\n";
        d_os << setw(7) << bit()[0] + offatom
                << setw(7) << bit()[1] + offatom
                << setw(5) << bit().type() + 1 << "\n";
        ++count;
      }
    }
    offatom += sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  // BONDANGLETYPE block
  num = gff.numAngleTypes();
  d_os << "BONDANGLETYPE\n"
          << "#  NTTY: number of bond angle types\n"
          << num << "\n"
          << "#  CT: force constant\n"
          << "#  T0: bond angle at minimum energy in degrees\n"
          << "#         CT          T0\n";

  for (int i = 0; i < num; ++i) {
    if (i > 0 && !(i % 10))d_os << "# " << i << "\n";
    d_os.precision(5);
    d_os.setf(ios::fixed, ios::floatfield);
    d_os.setf(ios::scientific, ios::floatfield);
    d_os << setw(12) << gff.angleType(i).fc()
            << setw(12) << gff.angleType(i).t0() << "\n";
  }
  d_os << "END\n";

  // BONDANGLEH & BONDANGLE block
  num = 0;
  for (int i = 0; i < sys.numMolecules(); ++i) {
    AngleIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
      if (sys.mol(i).topology().atom(bit()[0]).isH() ||
              sys.mol(i).topology().atom(bit()[1]).isH() ||
              sys.mol(i).topology().atom(bit()[2]).isH())
        ++num;
    }
  }
  d_os << "BONDANGLEH\n"
          << "#  NTHEH: number of bond angles involving H atoms in solute\n"
          << num << "\n"
          << "#  ITH, JTH, KTH: atom sequence numbers\n"
          << "#    of atoms forming a bond angle in solute\n"
          << "#  ICTH: bond angle type code\n"
          << "#   ITH    JTH    KTH ICTH\n";

  for (int i = 0, offatom = 1; i < sys.numMolecules(); ++i) {
    AngleIterator bit(sys.mol(i).topology());
    for (int count = 0; bit; ++bit) {
      if (sys.mol(i).topology().atom(bit()[0]).isH() ||
              sys.mol(i).topology().atom(bit()[1]).isH() ||
              sys.mol(i).topology().atom(bit()[2]).isH()) {

        if (count > 0 && !(count % 10))d_os << "# " << count << "\n";
        d_os << setw(7) << bit()[0] + offatom
                << setw(7) << bit()[1] + offatom
                << setw(7) << bit()[2] + offatom
                << setw(5) << bit().type() + 1 << "\n";
        ++count;
      }
    }
    offatom += sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  num = 0;
  for (int i = 0; i < sys.numMolecules(); ++i) {
    AngleIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
      if (!sys.mol(i).topology().atom(bit()[0]).isH() &&
              !sys.mol(i).topology().atom(bit()[1]).isH() &&
              !sys.mol(i).topology().atom(bit()[2]).isH())
        ++num;
    }
  }
  d_os << "BONDANGLE\n"
          << "#  NTHE: number of bond angles NOT\n"
          << "#        involving H atoms in solute\n"
          << num << "\n"
          << "#  IT, JT, KT: atom sequence numbers of atoms\n"
          << "#     forming a bond angle\n"
          << "#  ICT: bond angle type code\n"
          << "#    IT     JT     KT  ICT\n";

  for (int i = 0, offatom = 1; i < sys.numMolecules(); ++i) {
    AngleIterator bit(sys.mol(i).topology());
    for (int count = 0; bit; ++bit) {
      if (!sys.mol(i).topology().atom(bit()[0]).isH() &&
              !sys.mol(i).topology().atom(bit()[1]).isH() &&
              !sys.mol(i).topology().atom(bit()[2]).isH()) {
        if (count > 0 && !(count % 10))d_os << "# " << count << "\n";
        d_os << setw(7) << bit()[0] + offatom
                << setw(7) << bit()[1] + offatom
                << setw(7) << bit()[2] + offatom
                << setw(5) << bit().type() + 1 << "\n";
        ++count;
      }
    }
    offatom += sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  // IMPDIHEDRALTYPE block
  num = gff.numImproperTypes();
  d_os << "IMPDIHEDRALTYPE\n"
          << "#  NQTY: number of improper dihedrals\n"
          << num << "\n"
          << "#  CQ: force constant of improper dihedral per degrees square\n"
          << "#  Q0: improper dihedral angle at minimum energy in degrees\n"
          << "#         CQ          Q0\n";

  for (int i = 0; i < num; ++i) {
    if (i > 0 && !(i % 10))d_os << "# " << i << "\n";
    d_os.precision(5);
    d_os.setf(ios::fixed, ios::floatfield);
    d_os.setf(ios::scientific, ios::floatfield);
    d_os << setw(12) << gff.improperType(i).fc()
            << setw(12) << gff.improperType(i).q0() << "\n";
  }
  d_os << "END\n";

  // IMPDIHEDRALH & IMPDIHEDRAL block
  num = 0;
  for (int i = 0; i < sys.numMolecules(); ++i) {
    ImproperIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
      if (sys.mol(i).topology().atom(bit()[0]).isH() ||
              sys.mol(i).topology().atom(bit()[1]).isH() ||
              sys.mol(i).topology().atom(bit()[2]).isH() ||
              sys.mol(i).topology().atom(bit()[3]).isH())
        ++num;
    }
  }
  d_os << "IMPDIHEDRALH\n"
          << "#  NQHIH: number of improper dihedrals\n"
          << "#         involving H atoms in the solute\n"
          << num << "\n"
          << "#  IQH,JQH,KQH,LQH: atom sequence numbers\n"
          << "#     of atoms forming an improper dihedral\n"
          << "#  ICQH: improper dihedral type code\n"
          << "#   IQH    JQH    KQH    LQH ICQH\n";

  for (int i = 0, offatom = 1; i < sys.numMolecules(); ++i) {
    ImproperIterator bit(sys.mol(i).topology());
    for (int count = 0; bit; ++bit) {
      if (sys.mol(i).topology().atom(bit()[0]).isH() ||
              sys.mol(i).topology().atom(bit()[1]).isH() ||
              sys.mol(i).topology().atom(bit()[2]).isH() ||
              sys.mol(i).topology().atom(bit()[3]).isH()) {
        if (count > 0 && !(count % 10))d_os << "# " << count << "\n";
        d_os << setw(7) << bit()[0] + offatom
                << setw(7) << bit()[1] + offatom
                << setw(7) << bit()[2] + offatom
                << setw(7) << bit()[3] + offatom
                << setw(5) << bit().type() + 1 << "\n";
        ++count;
      }
    }
    offatom += sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  num = 0;
  for (int i = 0; i < sys.numMolecules(); ++i) {
    ImproperIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
      if (!sys.mol(i).topology().atom(bit()[0]).isH() &&
              !sys.mol(i).topology().atom(bit()[1]).isH() &&
              !sys.mol(i).topology().atom(bit()[2]).isH() &&
              !sys.mol(i).topology().atom(bit()[3]).isH())
        ++num;
    }
  }
  d_os << "IMPDIHEDRAL\n"
          << "#  NQHI: number of improper dihedrals NOT\n"
          << "#    involving H atoms in solute\n"
          << num << "\n"
          << "#  IQ,JQ,KQ,LQ: atom sequence numbers of atoms\n"
          << "#    forming an improper dihedral\n"
          << "#  ICQ: improper dihedral type code\n"
          << "#    IQ     JQ     KQ     LQ  ICQ\n";

  for (int i = 0, offatom = 1; i < sys.numMolecules(); ++i) {
    ImproperIterator bit(sys.mol(i).topology());
    for (int count = 0; bit; ++bit) {
      if (!sys.mol(i).topology().atom(bit()[0]).isH() &&
              !sys.mol(i).topology().atom(bit()[1]).isH() &&
              !sys.mol(i).topology().atom(bit()[2]).isH() &&
              !sys.mol(i).topology().atom(bit()[3]).isH()) {
        if (count > 0 && !(count % 10))d_os << "# " << count << "\n";
        d_os << setw(7) << bit()[0] + offatom
                << setw(7) << bit()[1] + offatom
                << setw(7) << bit()[2] + offatom
                << setw(7) << bit()[3] + offatom
                << setw(5) << bit().type() + 1 << "\n";
        ++count;
      }
    }
    offatom += sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  // DIHEDRALTYPE block
  num = gff.numDihedralTypes();

  d_os << "DIHEDRALTYPE\n"
          << "#  NPTY: number of dihedral types\n"
          << num << "\n"
          << "#  CP: force constant\n"
          << "#  PD: cosine of the phase shift\n"
          << "#  NP: multiplicity\n"
          << "#       CP        PD  NP\n";

  for (int i = 0; i < num; ++i) {
    if (i > 0 && !(i % 10))d_os << "# " << i << "\n";
    d_os.precision(5);
    d_os.setf(ios::fixed, ios::floatfield);
    d_os << setw(10) << gff.dihedralType(i).fc()
            << setw(10) << gff.dihedralType(i).pd()
            << setw(4) << gff.dihedralType(i).np() << "\n";
  }
  d_os << "END\n";

  // DIHEDRALH & DIHEDRAL block
  num = 0;
  for (int i = 0; i < sys.numMolecules(); ++i) {
    DihedralIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
      if (sys.mol(i).topology().atom(bit()[0]).isH() ||
              sys.mol(i).topology().atom(bit()[1]).isH() ||
              sys.mol(i).topology().atom(bit()[2]).isH() ||
              sys.mol(i).topology().atom(bit()[3]).isH())
        ++num;
    }
  }
  d_os << "DIHEDRALH\n"
          << "#  NPHIH: number of dihedrals involving H atoms in solute\n"
          << num << "\n"
          << "#  IPH, JPH, KPH, LPH: atom sequence numbers\n"
          << "#    of atoms forming a dihedral\n"
          << "#  ICPH: dihedral type code\n"
          << "#   IPH    JPH    KPH    LPH ICPH\n";

  for (int i = 0, offatom = 1; i < sys.numMolecules(); ++i) {
    DihedralIterator bit(sys.mol(i).topology());
    for (int count = 0; bit; ++bit) {
      if (sys.mol(i).topology().atom(bit()[0]).isH() ||
              sys.mol(i).topology().atom(bit()[1]).isH() ||
              sys.mol(i).topology().atom(bit()[2]).isH() ||
              sys.mol(i).topology().atom(bit()[3]).isH()) {
        if (count > 0 && !(count % 10))d_os << "# " << count << "\n";
        d_os << setw(7) << bit()[0] + offatom
                << setw(7) << bit()[1] + offatom
                << setw(7) << bit()[2] + offatom
                << setw(7) << bit()[3] + offatom
                << setw(5) << bit().type() + 1 << "\n";
        ++count;
      }
    }
    offatom += sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  num = 0;
  for (int i = 0; i < sys.numMolecules(); ++i) {
    DihedralIterator bit(sys.mol(i).topology());
    for (; bit; ++bit) {
      if (!sys.mol(i).topology().atom(bit()[0]).isH() &&
              !sys.mol(i).topology().atom(bit()[1]).isH() &&
              !sys.mol(i).topology().atom(bit()[2]).isH() &&
              !sys.mol(i).topology().atom(bit()[3]).isH())
        ++num;
    }
  }
  d_os << "DIHEDRAL\n"
          << "#  NPHI: number of dihedrals NOT involving H atoms in solute\n"
          << num << "\n"
          << "#  IP, JP, KP, LP: atom sequence numbers\n"
          << "#     of atoms forming a dihedral\n"
          << "#  ICP: dihedral type code\n"
          << "#    IP     JP     KP     LP  ICP\n";

  for (int i = 0, offatom = 1; i < sys.numMolecules(); ++i) {
    DihedralIterator bit(sys.mol(i).topology());
    for (int count = 0; bit; ++bit) {
      if (!sys.mol(i).topology().atom(bit()[0]).isH() &&
              !sys.mol(i).topology().atom(bit()[1]).isH() &&
              !sys.mol(i).topology().atom(bit()[2]).isH() &&
              !sys.mol(i).topology().atom(bit()[3]).isH()) {
        if (count > 0 && !(count % 10))d_os << "# " << count << "\n";
        d_os << setw(7) << bit()[0] + offatom
                << setw(7) << bit()[1] + offatom
                << setw(7) << bit()[2] + offatom
                << setw(7) << bit()[3] + offatom
                << setw(5) << bit().type() + 1 << "\n";
        ++count;
      }
    }
    offatom += sys.mol(i).numAtoms();
  }
  d_os << "END\n";

  // LJPARAMETERS block
  num = gff.numLJTypes();
  d_os << "LJPARAMETERS\n"
          << "#  NRATT2: number of LJ interaction types = NRATT*(NRATT+1)/2\n"
          << num << "\n"
          << "#  IAC,JAC: integer (van der Waals) atom type code\n"
          << "#  C12: r**(-12) term in nonbonded interactions\n"
          << "#   C6: r**(-6) term in nonbonded interactions\n"
          << "# CS12: r**(-12) term in 1-4 nonbonded interactions\n"
          << "#  CS6: r**(-6) term in 1-4 nonbonded interactions\n"
          << "# IAC  JAC           C12            C6          CS12           CS6\n";

  for (int i = 0; i < gff.numAtomTypeNames(); ++i) {
    for (int j = 0; j <= i; ++j) {
      d_os.precision(6);
      d_os.setf(ios::fixed, ios::floatfield);
      d_os.setf(ios::scientific, ios::floatfield);
      LJType lj(gff.ljType(AtomPair(i, j)));
      d_os << setw(5) << j + 1
              << setw(5) << i + 1
              << setw(14) << lj.c12()
              << setw(14) << lj.c6()
              << setw(14) << lj.cs12()
              << setw(14) << lj.cs6() << "\n";
    }
    d_os << "#\n";
  }
  d_os << "END\n";

  //SOLVENTATOM block
  d_os << "SOLVENTATOM\n"
          << "#  NRAM: number of atoms per solvent molecule\n"
          << sys.sol(0).topology().numAtoms() << "\n"
          << "#     I: solvent atom sequence number\n"
          << "#  IACS: integer (van der Waals) atom type code\n"
          << "#  ANMS: atom name of solvent atom\n"
          << "#  MASS: mass of solvent atom\n"
          << "#   CGS: charge of solvent atom\n"
          << "#  I  ANMS IACS      MASS        CGS\n";

  for (int j = 0; j < sys.sol(0).topology().numAtoms(); ++j) {
    d_os.precision(5);
    d_os.setf(ios::fixed, ios::floatfield);
    d_os << setw(4) << j + 1 << " "
            << setw(5) << sys.sol(0).topology().atom(j).name()
            << setw(4) << sys.sol(0).topology().atom(j).iac() + 1
            << setw(11) << sys.sol(0).topology().atom(j).mass()
            << setw(11) << sys.sol(0).topology().atom(j).charge();
    d_os << "\n";
  }
  d_os << "END\n";

  //SOLVENTCONSTR bock
  num = 0;
  ConstraintIterator dit(sys.sol(0).topology());
  for (; dit; ++dit) ++num;
  d_os << "SOLVENTCONSTR\n"
          << "#  NCONS: number of constraints\n"
          << num << "\n"
          << "#  ICONS, JCONS: atom sequence numbers forming constraint\n"
          << "#   CONS constraint length\n"
          << "#ICONS JCONS         CONS\n";

  ConstraintIterator cit(sys.sol(0).topology());

  for (; cit; ++cit) {
    d_os.precision(7);
    d_os.setf(ios::fixed, ios::floatfield);

    d_os << setw(5) << cit()[0] + 1
            << setw(5) << cit()[1] + 1
            << setw(15) << cit().dist() << "\n";
  }
  d_os << "END\n";
  d_os << "# end of topology file" << endl;
}


