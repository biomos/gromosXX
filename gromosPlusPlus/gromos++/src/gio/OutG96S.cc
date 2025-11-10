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

// gio_OutG96S.cc
#include "OutG96S.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <string>

#include "../gromos/Exception.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/Solvent.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gmath/Vec.h"
#include "../gmath/Matrix.h"
#include "../gcore/Box.h"
#include "../gcore/Remd.h"
#include "../args/Arguments.h"
#include "../utils/AtomSpecifier.h"
#include "OutCoordinates.h"


using gio::OutG96S;
using namespace gcore;
using namespace std;
using namespace utils;

class gio::OutG96S_i {
  friend class gio::OutG96S;
  ostream &d_os;
  bool posres;
  int d_count, d_res_off, d_switch;

  OutG96S_i(ostream &os, bool posres = false) :
  d_os(os), posres(posres), d_count(0), d_res_off(1), d_switch() {
    d_switch = 0;
  }

  ~OutG96S_i() {
  }

  void select(const string &thing);
  void writeSingleM(const Molecule &mol);
  void writeSingleV(const gcore::System &sys);
  void writeSingleS(const Solvent &sol);
  void writeSingleM_vel(const Molecule &mol);
  void writeSingleS_vel(const Solvent &sol);
  void writeBox(const Box &box);
  void writeTriclinicBox(const Box &box);
  void writeGenBox(const Box &box);
  void writeRemd(const Remd &remd);
  void writeAtomSpecifier(const AtomSpecifier & atoms, bool vel = false);
};

OutG96S::OutG96S(ostream &os, bool p) :
OutCoordinates(),
d_this(new OutG96S_i(os, p)), posres(p) {
}

OutG96S::OutG96S(bool p) :
OutCoordinates(), posres(p) {
  d_this = 0;
}

OutG96S::~OutG96S() {
  if (d_this)delete d_this;
}

void OutG96S::writeTitle(const string &title) {
  d_this->d_os << "TITLE\n" << title << "\nEND\n";
}

void OutG96S::writeTimestep(const int step, const double time) {
  d_this->d_os.precision(9);
  d_this->d_os.setf(std::ios::fixed, std::ios::floatfield);

  d_this->d_os << "TIMESTEP\n"
          << std::setw(18)
          << step
          << std::setw(20)
          << time
          << "\n#if @time flag is used the value for step refers to the"
          << "\n#step-th configuration in the original trajectory file"          
          << "\nEND\n";
}

void OutG96S::select(const string &thing) {
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

void OutG96S::open(ostream &os) {
  if (d_this) {
    delete d_this;
  }
  d_this = new OutG96S_i(os, posres);
}

void OutG96S::close() {
  if (d_this)delete d_this;
  d_this = 0;
}

OutG96S &OutG96S::operator<<(const gcore::System &sys) {
  // what do we have to write out
  bool writePos = false;
  bool writeVel = false;
  for (int m = 0; m < sys.numMolecules(); m++) {
    if (sys.mol(m).numPos()) {
      writePos = true;
    }
    if (sys.mol(m).numVel()) {
      writeVel = true;
    }
  }
  for (int s = 0; s < sys.numSolvents(); s++) {
    if (sys.sol(s).numPos()) {
      writePos = true;
    }
    if (sys.sol(s).numVel()) {
      writeVel = true;
    }
  }

  if (!posres) {
    if (sys.hasRemd) {
      d_this->d_os << "REMD\n";
      d_this->writeRemd(sys.remd());
      d_this->d_os << "END\n";
    }
  }

  if (writePos) {
    d_this->d_count = 0;
    d_this->d_res_off = 1;
    if (d_this->d_switch == 2) {
      d_this->d_res_off = 0;
    }
    if (!posres) {
      d_this->d_os << "POSITION\n";
    } else {
      d_this->d_os << "POSRESSPEC\n";
    }
    if (d_this->d_switch == 0 || d_this->d_switch == 1 || 
        d_this->d_switch == 3 || d_this->d_switch == 4)
      for (int i = 0; i < sys.numMolecules(); ++i)
        d_this->writeSingleM(sys.mol(i));
    if (d_this->d_switch == 3 || d_this->d_switch == 4 || 
        d_this->d_switch == 5) {
      d_this->writeSingleV(sys);}
    if (d_this->d_switch == 1 || d_this->d_switch == 2 ||
        d_this->d_switch == 4 || d_this->d_switch == 5) {
      for (int i = 0; i < sys.numSolvents(); ++i)
        d_this->writeSingleS(sys.sol(i));
    }
    d_this->d_os << "END\n";
  }

  if (!posres) {
    if (writeVel) {
      d_this->d_count = 0;
      d_this->d_res_off = 1;
      if (d_this->d_switch == 2) {
        d_this->d_res_off = 0;
      }
      d_this->d_os << "VELOCITY\n";
      if (d_this->d_switch <= 1) {
        for (int i = 0; i < sys.numMolecules(); ++i) {
          d_this->writeSingleM_vel(sys.mol(i));
        }
      }
      if (d_this->d_switch >= 1) {
        for (int i = 0; i < sys.numSolvents(); ++i) {
          d_this->writeSingleS_vel(sys.sol(i));
        }
      }
      d_this->d_os << "END\n";
    }

    if (args::Arguments::outG96) {
      switch (sys.box().boxformat()) {
        case gcore::Box::box96:
        {
          d_this->d_os << "BOX\n";
          d_this->writeBox(sys.box());
          d_this->d_os << "END\n";
          break;
        }
        case gcore::Box::triclinicbox:
        {
          d_this->d_os << "TRICLINICBOX\n";
          d_this->writeTriclinicBox(sys.box());
          d_this->d_os << "END\n";
          break;
        }
        case gcore::Box::genbox:
        {
          d_this->d_os << "GENBOX\n";
          d_this->writeGenBox(sys.box());
          d_this->d_os << "END\n";
          break;
        }
        default:
          throw gromos::Exception("OutG96", "Don't know how to handle boxformat");
      }
    } else {
      d_this->d_os << "GENBOX\n";
      d_this->writeGenBox(sys.box());
      d_this->d_os << "END\n";
    }
  }

  return *this;
}

OutG96S &OutG96S::operator<<(const AtomSpecifier & atoms) {
  const System & sys = *(atoms.sys());
  // what do we have to write out
  bool writePos = false;
  bool writeVel = false;
  for (int m = 0; m < sys.numMolecules(); m++) {
    if (sys.mol(m).numPos()) {
      writePos = true;
    }
    if (sys.mol(m).numVel()) {
      writeVel = true;
    }
  }
  for (int s = 0; s < sys.numSolvents(); s++) {
    if (sys.sol(s).numPos()) {
      writePos = true;
    }
    if (sys.sol(s).numVel()) {
      writeVel = true;
    }
  }

  if (!posres) {
    if (sys.hasRemd) {
      d_this->d_os << "REMD\n";
      d_this->writeRemd(sys.remd());
      d_this->d_os << "END\n";
    }
  }

  if (writePos) {
    d_this->writeAtomSpecifier(atoms);
  }

  if (!posres) {
    if (writeVel) {
      d_this->writeAtomSpecifier(atoms, true);
    }

    if (args::Arguments::outG96) {
      switch (sys.box().boxformat()) {
        case gcore::Box::box96:
        {
          d_this->d_os << "BOX\n";
          d_this->writeBox(sys.box());
          d_this->d_os << "END\n";
          break;
        }
        case gcore::Box::triclinicbox:
        {
          d_this->d_os << "TRICLINICBOX\n";
          d_this->writeTriclinicBox(sys.box());
          d_this->d_os << "END\n";
          break;
        }
        case gcore::Box::genbox:
        {
          d_this->d_os << "GENBOX\n";
          d_this->writeGenBox(sys.box());
          d_this->d_os << "END\n";
          break;
        }
        default:
          throw gromos::Exception("OutG96", "Don't know how to handle boxformat");
      }
    } else {
      d_this->d_os << "GENBOX\n";
      d_this->writeGenBox(sys.box());
      d_this->d_os << "END\n";
    }
  }

  return *this;
}

void gio::OutG96S_i::writeRemd(const Remd &remd) {
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.precision(6);
  d_os << setw(15) << remd.id()
          << setw(10) << remd.run()
          << setprecision(1) << setw(10) << remd.temperature()
          << setprecision(6) << setw(10) << remd.lambda()
          << "\n"
          << setw(15) << remd.Ti()
          << setw(10) << remd.li()
          << setw(10) << remd.Tj()
          << setw(10) << remd.lj()
          << setw(10) << remd.reeval()
          << "\n";
}

void gio::OutG96S_i::writeBox(const Box &box) {
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.precision(9);

  d_os << setw(15) << box.K()[0]
          << setw(15) << box.L()[1]
          << setw(15) << box.M()[2] << endl;
}

void gio::OutG96S_i::writeTriclinicBox(const Box &box) {
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.precision(9);

  d_os << setw(8) << box.ntb() << endl;
  for (int i = 0; i < 3; ++i) {
    d_os << setw(15) << box.K()[i]
            << setw(15) << box.L()[i]
            << setw(15) << box.M()[i] << endl;
  }

}

void gio::OutG96S_i::writeSingleM(const Molecule &mol) {
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.setf(ios::unitbuf);
  d_os.precision(9);
  for (int i = 0; i < mol.numAtoms(); ++i) {
    ++d_count;
    int res = mol.topology().resNum(i);
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(5) << res + d_res_off;
    d_os.setf(ios::left, ios::adjustfield);
    d_os << ' ' << setw(6) << mol.topology().resName(res).c_str()
            << setw(6) << mol.topology().atom(i).name().c_str();
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(6) << d_count;
    if (posres) {
      d_os << endl;
    } else {
      d_os << setw(15) << mol.pos(i)[0]
              << setw(15) << mol.pos(i)[1]
              << setw(15) << mol.pos(i)[2] << endl;
    }
  }
  d_res_off += mol.topology().numRes();
}

void gio::OutG96S_i::writeSingleV(const gcore::System &sys) {
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.precision(9);

  int res=0;
  int mn=1;
  if(d_count!=0){ 
    mn=sys.numMolecules();
  }

  for (int i = 0; i < sys.vas().numVirtualAtoms(); ++i) {
    ++d_count;
    //cout << "Virtual";
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(5) << res + d_res_off;
    d_os.setf(ios::left, ios::adjustfield);
    d_os << ' ' << setw(6) << "VIR"
            << "VIRT";
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(8) << d_count;
    if (posres) {
      d_os << endl;
    } else {
      d_os << setw(15) << sys.vas().atom(i).pos()[0]
              << setw(15) << sys.vas().atom(i).pos()[1]
              << setw(15) << sys.vas().atom(i).pos()[2] << endl;
    }
  }
}

void gio::OutG96S_i::writeSingleS(const Solvent &sol) {

  int na = sol.topology().numAtoms();
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.setf(ios::unitbuf);
  d_os.precision(9);
  for (int i = 0; i < sol.numPos(); ++i) {
    ++d_count;
    int res = i / na;
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(5) << res + 1; // /+ d_res_off;
    d_os.setf(ios::left, ios::adjustfield);
    d_os << ' ' << setw(6) << sol.topology().solvName().c_str()
            << setw(6) << sol.topology().atom(i % na).name().c_str();
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(6) << d_count;
    if (posres) {
      d_os << endl;
    } else {
      d_os << setw(15) << sol.pos(i)[0]
              << setw(15) << sol.pos(i)[1]
              << setw(15) << sol.pos(i)[2] << endl;
    }
  }
  //d_res_off += sol.numPos()/na;
}

void gio::OutG96S_i::writeSingleM_vel(const Molecule &mol) {
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.setf(ios::unitbuf);
  d_os.precision(9);
  for (int i = 0; i < mol.numVel(); ++i) {
    ++d_count;
    int res = mol.topology().resNum(i);
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(5) << res + d_res_off;
    d_os.setf(ios::left, ios::adjustfield);
    d_os << ' ' << setw(6) << mol.topology().resName(res).c_str()
            << setw(6) << mol.topology().atom(i).name().c_str();
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(6) << d_count
            << setw(15) << mol.vel(i)[0]
            << setw(15) << mol.vel(i)[1]
            << setw(15) << mol.vel(i)[2] << endl;
  }
  d_res_off += mol.topology().numRes();
}

void gio::OutG96S_i::writeSingleS_vel(const Solvent &sol) {
  int na = sol.topology().numAtoms();
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.setf(ios::unitbuf);
  d_os.precision(9);
  for (int i = 0; i < sol.numVel(); ++i) {
    ++d_count;
    int res = i / na;
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(5) << res + 1; // + d_res_off;
    d_os.setf(ios::left, ios::adjustfield);
    d_os << ' ' << setw(6) << sol.topology().solvName().c_str()
            << setw(6) << sol.topology().atom(i % na).name().c_str();
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(6) << d_count
            << setw(15) << sol.vel(i)[0]
            << setw(15) << sol.vel(i)[1]
            << setw(15) << sol.vel(i)[2] << endl;
  }
  //d_res_off += sol.numVel()/na;
}

void gio::OutG96S_i::writeGenBox(const Box &box) {
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.precision(9);
  const double k = box.K().abs();
  const double l = box.L().abs();
  const double m = box.M().abs();
  d_os << setw(8) << box.ntb() << endl;
  if (box.ntb() == gcore::Box::vacuum) {
    d_os << setw(15) << 0.0 << setw(15) << 0.0 << setw(15) << 0.0 << endl
            << setw(15) << 0.0 << setw(15) << 0.0 << setw(15) << 0.0 << endl
            << setw(15) << 0.0 << setw(15) << 0.0 << setw(15) << 0.0 << endl
            << setw(15) << box.X() << setw(15) << box.Y() << setw(15) << box.Z() << endl;
  } else {
    d_os << setw(15) << k
            << setw(15) << l
            << setw(15) << m << endl;
    d_os << setw(15) << acos(box.L().dot(box.M()) / (l * m))*180 / M_PI
            << setw(15) << acos(box.K().dot(box.M()) / (k * m))*180 / M_PI
            << setw(15) << acos(box.K().dot(box.L()) / (k * l))*180 / M_PI << endl;

    // calculate the Euler rotation angles as described in Phils manuscript:
    // "GROMOS01: Description of the changes", Philippe Huenenberger, October 5, 2004
    gmath::Vec x = box.K().normalize();
    gmath::Vec y = (box.L() - (box.L().dot(x) * x)).normalize();
    gmath::Vec z = x.cross(y);

    gmath::Matrix R_(x, y, z);
    double R11R21 = R_(0,0) * R_(0,0) + R_(1,0) * R_(1,0);
    double theta, psi, phi;
    if(R11R21 == 0.0) {
        int sign = 1;
        if(R_(2,0)<0) sign = -1;
        theta = -sign*M_PI/2;
        psi = 0.0;
        sign = 1;
        if(R_(0,1)<0) sign = -1;
        phi = -sign*acos(R_(1,1));
    } else {
        int sign =1;
        if(R_(2,0)<0) sign = -1;
        theta = -sign*acos(sqrt(R_(0,0)*R_(0,0)+R_(1,0)*R_(1,0)));
        sign = 1;
        if((R_(2,1)/cos(theta))<0) sign = -1;
        psi = sign*acos(R_(2,2)/cos(theta));
        sign = 1;
        if((R_(1,0)/cos(theta))<0) sign = -1;
        phi = sign*acos(R_(0,0)/cos(theta));
    }
    
    d_os << setw(15) << phi/M_PI*180
	 << setw(15) << theta/M_PI*180
	 << setw(15) << psi/M_PI*180 << endl;

    d_os << setw(15) << box.X()
            << setw(15) << box.Y()
            << setw(15) << box.Z() << endl;
  }
}

void gio::OutG96S_i::writeAtomSpecifier(const AtomSpecifier & atoms, bool vel) {
  if (posres)
    d_os << "POSRESSPEC" << endl;
  else if (vel)
    d_os << "VELOCITY" << endl;
  else
    d_os << "POSITION" << endl;

  d_os.setf(ios::fixed, ios::floatfield);
  d_os.precision(9);
  d_os << "# selected " << atoms.size() << " atoms" << endl;
  d_os.setf(ios::unitbuf);

  for (unsigned int i = 0; i < atoms.size(); ++i) {
    d_os.setf(ios::right, ios::adjustfield);
    int offset = 1;
    if (atoms.mol(i) >= 0) {
      for (int j = 0; j < atoms.mol(i); ++j)
        offset += atoms.sys()->mol(j).topology().numRes();
    }
    d_os << setw(5) << atoms.resnum(i) + offset;
    d_os.setf(ios::left, ios::adjustfield);
    string res = atoms.resname(i);

    if (atoms.mol(i) < 0) res = "SOLV";
    d_os << ' ' << setw(6) << res
            << setw(6) << atoms.name(i);
    d_os.setf(ios::right, ios::adjustfield);
    d_os << setw(6) << atoms.gromosAtom(i) + 1;
    if (posres) {
      d_os << endl;
    } else {
      if (!vel) {
        d_os << setw(15) << atoms.pos(i)[0]
                << setw(15) << atoms.pos(i)[1]
                << setw(15) << atoms.pos(i)[2] << endl;
      } else {
        d_os << setw(15) << atoms.vel(i)[0]
                << setw(15) << atoms.vel(i)[1]
                << setw(15) << atoms.vel(i)[2] << endl;
      }
    }
  }
  d_os << "END" << endl;
}

