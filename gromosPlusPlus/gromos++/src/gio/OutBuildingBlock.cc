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
#include "OutBuildingBlock.h"

#include <cassert>
#include <iostream>
#include <iomanip>
#include <set>
#include <string>

#include "../gcore/BuildingBlock.h"
#include "../gcore/LJException.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Exclusion.h"
#include "../gcore/BbSolute.h"
#include "../gcore/Bond.h"
#include "../gcore/Angle.h"
#include "../gcore/Dihedral.h"
#include "../gcore/Improper.h"
#include "../gcore/Constraint.h"
#include "../gcore/Exclusion.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/SolventTopology.h"

using namespace std;
using namespace gcore;
using namespace gio;

gio::OutBuildingBlock::OutBuildingBlock(std::ostream &os) : 
        d_title("No title"), d_os(os) {
}

gio::OutBuildingBlock::~OutBuildingBlock() {
}

void gio::OutBuildingBlock::setTitle(const std::string& title) {
  d_title = title;
}

void gio::OutBuildingBlock::write(const gcore::BuildingBlock& bb) {
  d_os << "TITLE" << endl
          << d_title << endl
          << "END" << endl;
  d_os.precision(8);
  d_os << "PHYSICALCONSTANTS\n"
          << "# FPEPSI: 1.0/(4.0*PI*EPS0) (EPS0 is the permittivity of vacuum)\n"
          << setw(15) << bb.Fpepsi()
          << "\n# HBAR: Planck's constant HBAR = H/(2* PI)\n"
          << setw(15) << bb.Hbar()
          << "\n# SPDL: Speed of light (nm/ps)\n"
          << setw(15) << bb.Spdl()
          << "\n# BOLTZ: Boltzmann's constant kB\n"
          << setw(15) << bb.Boltz()
          << "\nEND\n";
  d_os << "LINKEXCLUSIONS" << endl
          << "# nearest neighbour exclusions when linking" << endl
          << "# NRNE" << endl
          << setw(5) << bb.LinkExclusions() << endl
          << "END" << endl;

  for (int i = 0; i < bb.numBbSolutes(); ++i) {
    writeSingle(bb.bb(i), BBTypeSolute);
  }
  for (int i = 0; i < bb.numBbSolvents(); ++i) {
    writeSolvent(bb.bs(i));
  }
  for (int i = 0; i < bb.numBbEnds(); ++i) {
    writeSingle(bb.be(i), BBTypeEnd);
  }
}

void gio::OutBuildingBlock::writeSingle(const gcore::BbSolute& bb, BBType type) {
  int last_few = 0;
  bool endgroup = type == BBTypeEnd;

  if (endgroup) {
    d_os << "MTBUILDBLEND" << endl;
    last_few = bb.rep();
  } else {
    d_os << "MTBUILDBLSOLUTE" << endl;
    last_few = bb.numPexcl();
  }
  d_os << "# building block (residue, nucleotide, etc.)" << endl;
  d_os << "# RNME" << endl;
  d_os << bb.resName() << endl;
  if (endgroup) {
    d_os << "# number of atoms, number of atoms to be replaced" << endl;
    d_os << "# NMAT,NREP" << endl;
  } else {
    d_os << "# number of atoms, number of preceding exclusions" << endl;
    d_os << "# NMAT,NLIN" << endl;
  }
  d_os << setw(5) << bb.numAtoms();
  if (endgroup)
    d_os << setw(5) << bb.rep() << endl;
  else {
    d_os << setw(5) << bb.numPexcl() << endl;
    d_os << "# preceding exclusions" << endl;
    d_os << "#ATOM                               MAE MSAE" << endl;
    for (int i = 0; i < bb.numPexcl(); i++) {
      d_os << setw(5) << i + 1 - bb.numPexcl()
              << setw(34) << bb.pexcl(i).size();
      for (int j = 0; j < bb.pexcl(i).size(); j++)
        d_os << setw(5) << bb.pexcl(i).atom(j) + 1;
      d_os << endl;
    }
  }

  d_os << "# atoms" << endl;
  d_os << "#ATOM ANM  IACM MASS        CGMICGM MAE MSAE" << endl;
  for (int i = 0; i < bb.numAtoms(); i++) {
    if (i == bb.numAtoms() - last_few) {
      if (endgroup)
        d_os << "# replacing atoms" << endl;
      else
        d_os << "# trailing atoms" << endl
              << "#ATOM ANM  IACM MASS        CGMICGM" << endl;
    }
    d_os << setw(5) << i + 1 << ' ';

    d_os.setf(ios::left, ios::adjustfield);

    d_os << setw(4) << bb.atom(i).name();
    d_os.setf(ios::fixed, ios::adjustfield);
    d_os.precision(5);
    d_os.setf(ios::fixed, ios::floatfield);

    d_os << setw(5) << bb.atom(i).iac() + 1
            << setw(5) << int(bb.atom(i).mass()) + 1
            << setw(11) << bb.atom(i).charge()
            << setw(4) << bb.atom(i).chargeGroup();

    if (i < bb.numAtoms() - last_few) {
      d_os << setw(4) << bb.atom(i).exclusion().size();
      for (int j = 0; j < bb.atom(i).exclusion().size(); j++) {
        d_os << setw(5) << bb.atom(i).exclusion().atom(j) + 1;
        if ((j + 1) % 6 == 0 && j + 1 < bb.atom(i).exclusion().size())
          d_os << endl << setw(39) << " ";
      }
    }
    d_os << endl;
  }
  d_os << "# bonds" << endl;
  d_os << "#  NB" << endl;
  int numBonds = 0;
  {
    BondIterator bi(bb);
    for (; bi; ++bi) numBonds++;
  }
  d_os << setw(5) << numBonds << endl;
  d_os << "#  IB   JB  MCB" << endl;
  BondIterator bi(bb);
  for (; bi; ++bi) {
    d_os << setw(5) << bi()[0] + 1
            << setw(5) << bi()[1] + 1
            << setw(5) << bi().type() + 1 << endl;
  }
  d_os << "# bond angles" << endl;
  d_os << "# NBA" << endl;
  int numAngles = 0;
  {
    AngleIterator ai(bb);
    for (; ai; ++ai) numAngles++;
  }
  d_os << setw(5) << numAngles << endl;
  d_os << "#  IB   JB   KB  MCB" << endl;

  AngleIterator ai(bb);
  for (; ai; ++ai) {
    d_os << setw(5) << ai()[0] + 1
            << setw(5) << ai()[1] + 1
            << setw(5) << ai()[2] + 1
            << setw(5) << ai().type() + 1 << endl;
  }
  d_os << "# improper dihedrals" << endl;
  d_os << "# NIDA" << endl;
  int numImpropers = 0;
  {
    ImproperIterator ii(bb);
    for (; ii; ++ii) numImpropers++;
  }

  d_os << setw(5) << numImpropers << endl;
  d_os << "#  IB   JB   KB   LB  MCB" << endl;
  ImproperIterator ii(bb);
  for (; ii; ++ii) {
    d_os << setw(5) << ii()[0] + 1
            << setw(5) << ii()[1] + 1
            << setw(5) << ii()[2] + 1
            << setw(5) << ii()[3] + 1
            << setw(5) << ii().type() + 1 << endl;
  }
  d_os << "# dihedrals" << endl;
  d_os << "# NDA" << endl;
  int numDihedrals = 0;
  {
    DihedralIterator di(bb);
    for (; di; ++di) numDihedrals++;
  }

  d_os << setw(5) << numDihedrals << endl;
  d_os << "#  IB   JB   KB   LB  MCB" << endl;
  DihedralIterator di(bb);
  for (; di; ++di) {
    d_os << setw(5) << di()[0] + 1
            << setw(5) << di()[1] + 1
            << setw(5) << di()[2] + 1
            << setw(5) << di()[3] + 1
            << setw(5) << di().type() + 1 << endl;
  }
  d_os << "# LJ exceptions" << endl;
  d_os << "# NEX" << endl;
  int numLJExceptions = 0;
  {
    LJExceptionIterator li(bb);
    for(; li; ++li) numLJExceptions++;
  }
  d_os << setw(5) << numLJExceptions << endl;
  d_os << "#  IB   JB  MCB  NCO  IND  CON" << endl;
  LJExceptionIterator li(bb);
  for(; li; ++li) {
    d_os << setw(5) << li()[0] + 1
            << setw(5) << li()[1] + 1
            << setw(5) << li().type() + 1
            << setw(5) << li().numcond();
    if (li().numcond()) {
      d_os << setw(5) << li().indicate();
      set<int>::const_iterator it = li().cond().begin(),
              to = li().cond().end();
      for(; it != to; ++it) {
        d_os << setw(5) << *it + 1;
      }
    }
    d_os << endl;
  }
  d_os << "END" << endl;
}

void gio::OutBuildingBlock::writeSolvent(const gcore::SolventTopology& bb) {
  d_os << "MTBUILDBLSOLVENT" << endl
          << "#solvent name" << endl
          << "#RNMES" << endl
          << bb.solvName() << endl
          << "#number of atoms" << endl
          << bb.numAtoms() << endl
          << "#atoms" << endl
          << "#ATOM ANM    IAC MASS    CG" << endl;
  for(int i = 0; i < bb.numAtoms(); ++i) {
    d_os << setw(5) << i + 1 << ' ';

    d_os.setf(ios::left, ios::adjustfield);

    d_os << setw(5) << bb.atom(i).name();
    d_os.setf(ios::fixed, ios::adjustfield);
    d_os.precision(5);
    d_os.setf(ios::fixed, ios::floatfield);

    d_os << setw(5) << bb.atom(i).iac() + 1
            << setw(5) << int(bb.atom(i).mass()) + 1
            << setw(11) << bb.atom(i).charge() << endl;
  }
  d_os << "# constraints" << endl;
  int numConstraints = 0;
  {
    ConstraintIterator ci(bb);
    for(; ci; ++ci) numConstraints++;
  }
  d_os << setw(5) << numConstraints << endl
          << "#  IB   JB  LENGTH" << endl;
  ConstraintIterator ci(bb);
  for(; ci; ++ci) {
    d_os.precision(7);
    d_os.setf(ios::fixed, ios::floatfield);
    d_os << setw(5) << ci()[0] + 1
            << setw(5) << ci()[1] + 1
            << setw(15) << ci().dist() << endl;
  }
}

