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

// gio_OutCif.cc
#include "OutCif.h"

#include <iomanip>
#include <iostream>
#include <string>

#include "OutCoordinates.h"
#include "../gcore/Box.h"
#include "../gcore/Molecule.h"
#include "../gcore/Solvent.h"
#include "../gcore/System.h"
#include "../gmath/Vec.h"
#include "../utils/AtomSpecifier.h"

using namespace gio;
using namespace gcore;
using namespace utils;

class gio::OutCif_i {
  friend class gio::OutCif;
  ostream &d_os;
  int d_count, d_res_off, d_switch;
  double d_factor;
  bool d_renumber;

  OutCif_i(ostream &os)
      : d_os{os}, d_count{0}, d_res_off{1}, d_switch{0}, d_factor{10.0} {}

  ~OutCif_i() {}

  void writeSingleK();
  void writeSingleM(const Molecule &mol, const int mn);
  void writeSingleV(const gcore::System &sys);
  void writeSingleS(const Solvent &sol);
  void writeSymmetry(string data_name);
  void writeCell(const Box &box, string data_name);
  void writeAtomSpecifier(const AtomSpecifier &atoms);
};

OutCif::OutCif(double factor, bool renumber)
    : OutCoordinates(), d_this{nullptr}, factor{factor}, renumber{renumber} {}

OutCif::OutCif(ostream &os, double factor, bool renumber)
    : OutCoordinates(), factor{factor}, renumber{renumber},
      d_this{new OutCif_i(os)} {
  d_this->d_factor = factor;
  d_this->d_renumber = renumber;
}

OutCif::~OutCif() {
  if (d_this)
    delete d_this;
}

void OutCif::select(const string &thing) {
  if (thing == "ALL") {
    d_this->d_switch = 1;
  } else if (thing == "SOLVENT") {
    d_this->d_switch = 2;
  } else if (thing == "SOLUTEV") {
    d_this->d_switch = 3;
  } else if (thing == "ALLV") {
    d_this->d_switch = 4;
  } else if (thing == "SOLVENTV") {
    d_this->d_switch = 5;
  } else {
    d_this->d_switch = 0;
  }
}

void OutCif::open(ostream &os) {
  if (d_this)
    delete d_this;

  d_this = new OutCif_i(os);
  d_this->d_factor = factor;
  d_this->d_renumber = renumber;
}

void OutCif::close() {
  if (d_this)
    delete d_this;

  d_this = nullptr;
}

void OutCif::writeTitle(const string &title) {
  data_name = title;
  d_this->d_os << "data_" << data_name << endl;
  d_this->d_os << "#" << endl;
}

// TODO: Only wrote as it is a virtual function
void OutCif::writeTimestep(const int step, const double time) {
  d_this->d_os.precision(9);
  d_this->d_os.setf(std::ios::fixed, std::ios::floatfield);

  d_this->d_os << "#REMARK   1  TIMESTEP\t" << std::setw(18) << step
               << std::setw(20) << time << "\n";
}

OutCif &OutCif::operator<<(const gcore::System &sys) {

  if (sys.hasBox) {
    d_this->writeCell(sys.box(), data_name);
    d_this->writeSymmetry(data_name);
  }

  d_this->d_count = 0;
  d_this->d_res_off = 1;
  d_this->writeSingleK();
  if (d_this->d_switch == 0 || d_this->d_switch == 1 || d_this->d_switch == 3 ||
      d_this->d_switch == 4) {
    for (auto i = 0; i != sys.numMolecules(); ++i) {
      d_this->writeSingleM(sys.mol(i), i + 1);
    }
  }
  if (d_this->d_switch == 3 || d_this->d_switch == 4 || d_this->d_switch == 5) {
    d_this->writeSingleV(sys);
  }
  if (d_this->d_switch == 1 || d_this->d_switch == 2 || d_this->d_switch == 4 ||
      d_this->d_switch == 5) {
    for (auto i = 0; i != sys.numSolvents(); ++i) {
      d_this->writeSingleS(sys.sol(i));
    }
  }
  d_this->d_os << "#" << endl;

  // TODO: write connect

  return *this;
}

OutCif &OutCif::operator<<(const utils::AtomSpecifier &atoms) {

  if (atoms.sys()->hasBox) {
    d_this->writeCell(atoms.sys()->box(), data_name);
    d_this->writeSymmetry(data_name);
  }

  d_this->writeAtomSpecifier(atoms);
  d_this->d_os << "#" << endl;

  return *this;
}

void gio::OutCif_i::writeCell(const Box &box, string data_name) {
  if (box.ntb() == Box::vacuum || box.ntb() == Box::truncoct) {
    return;
  }

  ++d_count;
  d_os.setf(ios::unitbuf);

  d_os << "_cell.entry_id\t" << data_name << endl;
  d_os << "_cell.length_a\t";
  d_os.precision(3);
  d_os << box.K().abs() * d_factor << endl;
  d_os << "_cell.length_b\t";
  d_os.precision(3);
  d_os << box.L().abs() * d_factor << endl;
  d_os << "_cell.length_c\t";
  d_os.precision(3);
  d_os << box.M().abs() * d_factor << endl;
  d_os << "_cell.angle_alpha\t";
  d_os.precision(2);
  d_os << box.alpha() << endl;
  d_os << "_cell.angle_beta\t";
  d_os.precision(2);
  d_os << box.beta() << endl;
  d_os << "_cell.angle_gamma\t";
  d_os.precision(2);
  d_os << box.gamma() << endl;
  d_os << "#" << endl;
}

void gio::OutCif_i::writeSymmetry(string data_name) {
  d_os << "_symmetry.entry_id\t" << data_name << endl;
  d_os << "_symmetry.space_group_name_H-M\t" << "'P 1'" << endl;
  d_os << "#" << endl;
}

void gio::OutCif_i::writeSingleK() {
  d_os << "loop_" << endl;
  d_os << "_atom_site.group_PDB" << endl;
  d_os << "_atom_site.id" << endl;
  d_os << "_atom_site.type_symbol" << endl;
  d_os << "_atom_site.label_atom_id" << endl;
  // d_os << "_atom_site.label_alt_id" << endl;
  d_os << "_atom_site.label_comp_id" << endl;
  d_os << "_atom_site.label_asym_id" << endl;
  // d_os << "_atom_site.label_entity_id" << endl;
  d_os << "_atom_site.label_seq_id" << endl;
  // d_os << "_atom_site.pdbx_PDB_ins_code" << endl;
  d_os << "_atom_site.Cartn_x" << endl;
  d_os << "_atom_site.Cartn_y" << endl;
  d_os << "_atom_site.Cartn_z" << endl;
  d_os << "_atom_site.occupancy" << endl;
  d_os << "_atom_site.B_iso_or_equiv" << endl;
  // d_os << "_atom_site.Cartn_x_esd" << endl;
  // d_os << "_atom_site.Cartn_y_esd" << endl;
  // d_os << "_atom_site.Cartn_z_esd" << endl;
  // d_os << "_atom_site.occupancy_esd" << endl;
  // d_os << "_atom_site.B_iso_or_equiv_esd" << endl;
  // d_os << "_atom_site.pdbx_formal_charge" << endl;
  // d_os << "_atom_site.auth_seq_id" << endl;
  // d_os << "_atom_site.auth_comp_id" << endl;
  // d_os << "_atom_site.auth_asym_id" << endl;
  // d_os << "_atom_site.auth_atom_id" << endl;
  // d_os << "_atom_site.pdbx_PDB_model_num" << endl;
}

void gio::OutCif_i::writeSingleM(const Molecule &mol, const int mn) {

  double bfac = 0;
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.setf(ios::unitbuf);
  d_os.precision(3);

  for (auto i = 0; i != mol.numAtoms(); ++i) {
    ++d_count;
    int res = mol.topology().resNum(i);

    d_os.setf(ios::left, ios::adjustfield);
    d_os << setw(6) << "ATOM";

    d_os.setf(ios::right, ios::adjustfield);
    d_os << ' ' << setw(6) << d_count;

    d_os.setf(ios::left, ios::adjustfield);
    d_os << ' ' << mol.topology().atom(i).name().substr(0, 1).c_str();
    if (mol.topology().atom(i).name().length() == 4) {
      d_os << ' ' << setw(5)
           << mol.topology().atom(i).name().substr(0, 4).c_str();
    } else {
      d_os << "  " << setw(4)
           << mol.topology().atom(i).name().substr(0, 3).c_str();
    }
    d_os << setw(4) << mol.topology().resName(res).substr(0, 3).c_str();
    char chain = ('A' + mn - 1);
    if (chain < 'A' || chain > 'Z')
      chain = 'Z';
    d_os << ' ' << setw(1) << chain;

    d_os.setf(ios::right, ios::adjustfield);
    d_os << ' ' << setw(5) << res + d_res_off << "    " << setw(8)
         << mol.pos(i)[0] * d_factor << setw(8) << mol.pos(i)[1] * d_factor
         << setw(8) << mol.pos(i)[2] * d_factor;
    d_os << "  1.00" << setw(6) << setprecision(2) << bfac << setprecision(3)
         << endl;
  }
  if (!d_renumber)
    d_res_off += mol.topology().numRes();
}

void gio::OutCif_i::writeSingleV(const gcore::System &sys) {

  double bfac = 0;
  int res = 0;
  int mn = 1;
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.setf(ios::unitbuf);
  d_os.precision(3);
  if (d_count != 0) {
    mn = sys.numMolecules();
  }

  for (auto i = 0; i != sys.vas().numVirtualAtoms(); ++i) {
    ++d_count;

    d_os.setf(ios::left, ios::adjustfield);
    d_os << setw(6) << "ATOM";

    d_os.setf(ios::right, ios::adjustfield);
    d_os << ' ' << setw(6) << d_count;

    d_os.setf(ios::left, ios::adjustfield);
    d_os << " X";
    d_os << ' ' << setw(5) << "VIRT";
    d_os << ' ' << setw(4) << "VRT";
    char chain = ('A' + mn - 1);
    if (chain < 'A' || chain > 'Z')
      chain = 'Z';
    d_os << ' ' << setw(1) << chain;

    d_os.setf(ios::right, ios::adjustfield);
    d_os << ' ' << setw(5) << res + d_res_off << "    " << setw(8)
         << sys.vas().atom(i).pos()[0] * d_factor << setw(8)
         << sys.vas().atom(i).pos()[1] * d_factor << setw(8)
         << sys.vas().atom(i).pos()[2] * d_factor;
    d_os << "  1.00" << setw(6) << setprecision(2) << bfac << setprecision(3)
         << endl;
  }
}

void gio::OutCif_i::writeSingleS(const Solvent &sol) {

  int na = sol.topology().numAtoms();
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.setf(ios::unitbuf);
  d_os.precision(3);

  for (auto i = 0; i != sol.numPos(); ++i) {
    ++d_count;
    int res = i / na;
    int nameid = i % na;

    d_os.setf(ios::left, ios::adjustfield);
    d_os << setw(6) << "HETATM";

    d_os.setf(ios::right, ios::adjustfield);
    d_os << ' ' << setw(6) << d_count;

    d_os.setf(ios::left, ios::adjustfield);
    d_os << ' ' << sol.topology().atom(nameid).name().substr(0, 1).c_str();
    if (sol.topology().atom(nameid).name().length() == 4) {
      d_os << ' ' << setw(5)
           << sol.topology().atom(nameid).name().substr(0, 4).c_str();
    } else {
      d_os << "  " << setw(4)
           << sol.topology().atom(nameid).name().substr(0, 3).c_str();
    }
    d_os << setw(4) << sol.topology().solvName().c_str();
    d_os << " .";

    d_os.setf(ios::right, ios::adjustfield);
    d_os << ' ' << setw(5) << res + 1 << "    " << setw(8)
         << sol.pos(i)[0] * d_factor << setw(8) << sol.pos(i)[1] * d_factor
         << setw(8) << sol.pos(i)[2] * d_factor;
    d_os << "  1.00  0.00" << endl;
  }
  d_res_off += sol.numPos() / na;
}

void gio::OutCif_i::writeAtomSpecifier(const AtomSpecifier &atoms) {

  d_os.setf(ios::fixed, ios::floatfield);
  d_os.setf(ios::unitbuf);
  d_os.precision(3);
  int res = 0;
  double bfac = 0;
  int count = 0;
  int resoff = 0;
  char chain = 'A';

  for (auto i = 0; i != atoms.size(); ++i) {
    int maxmol = atoms.mol(i);
    if (maxmol < 0)
      maxmol = atoms.sys()->numMolecules();
    count = atoms.atom(i);
    resoff = 0;
    for (auto j = 0; j != maxmol; ++j) {
      count += atoms.sys()->mol(j).numAtoms();
      resoff += atoms.sys()->mol(j).topology().numRes();
    }

    if (atoms.mol(i) < 0) {
      res = atoms.atom(i) / atoms.sys()->sol(0).topology().numAtoms();
    } else {
      d_os << setw(4)
           << atoms.sys()
                  ->mol(atoms.mol(i))
                  .topology()
                  .resName(res)
                  .substr(0, 4)
                  .c_str();
    }

    if (atoms.mol(i) < 0 || atoms.sys()->mol(atoms.mol(i)).numBfac() <= 0) {
      bfac = 0;
    } else {
      bfac = atoms.sys()->mol(atoms.mol(i)).bfac(atoms.atom(i));
    }

    d_os.setf(ios::left, ios::adjustfield);
    d_os << setw(6) << "ATOM";

    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(7) << count + 1;

    d_os.setf(ios::left, ios::adjustfield);
    d_os << "  " << setw(4) << atoms.name(i).substr(0, 3).c_str();
    if (atoms.mol(i) < 0) {
      d_os << setw(4) << "SOLV";
    } else {
      d_os << setw(4)
           << atoms.sys()
                  ->mol(atoms.mol(i))
                  .topology()
                  .resName(res)
                  .substr(0, 4)
                  .c_str();
    }
    if (chain < 'A' || chain > 'Z')
      chain = 'Z';
    d_os << setw(1) << chain;
    d_os.setf(ios::right, ios::adjustfield);

    int resn = res + 1;
    if (!d_renumber)
      resn += resoff;
    d_os << setw(4) << resn << "    " << setw(8) << atoms.pos(i)[0] * d_factor
         << setw(8) << atoms.pos(i)[1] * d_factor << setw(8)
         << atoms.pos(i)[2] * d_factor;
    d_os << "  1.00" << setw(6) << setprecision(2) << bfac << setprecision(3)
         << endl;
    if (atoms.mol(i) != atoms.mol(i + 1)) {
      chain += 1;
    }
  }
}
