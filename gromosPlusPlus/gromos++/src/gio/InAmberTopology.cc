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

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   AmberTopology.cpp
 * Author: bschroed
 *
 * Created on March 7, 2018, 3:19 PM
 */

#include "InAmberTopology.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "Ginstream.h"
#include "../gcore/Angle.h"
#include "../gcore/AngleType.h"
#include "../gcore/AtomPair.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/Bond.h"
#include "../gcore/BondType.h"
#include "../gcore/Dihedral.h"
#include "../gcore/DihedralType.h"
#include "../gcore/Exclusion.h"
#include "../gcore/ImproperType.h"
#include "../gcore/LJType.h"
#include "../gcore/LinearTopology.h"
#include "../gmath/Physics.h"
#include "../gromos/Exception.h"
#include "../utils/StringOps.h"

//#define DebugDihedH      //debbuing by bschroed
//#define Debug_DihedType      //debbuing by bschroed
//#define DebugDihedwoH      //debbuing by bschroed

using namespace std;
using namespace gcore;

namespace gio {

AmberTopology::AmberTopology(const AmberTopology& orig) {
}

AmberTopology::AmberTopology(string s)  : d_blocks() {
  string topo = readAmber(s);
  stringstream amberTop(topo);
  this->stream(amberTop, false); //give stream to class
  this->init(); // moved init in the right scope, than works
}

AmberTopology::~AmberTopology() {
  this->close();
}

string AmberTopology::readAmber(string inputfile) {
  try {
    ifstream amberTopFile(inputfile);
    stringstream newlines;
    bool atomname = false;
    bool starting = true;

    for (string line; std::getline(amberTopFile, line);) {
      //check for new blocks
      if (line.find("%FLAG") != string::npos) { // new block
        if (starting) {
          newlines << utils::StringOps::replaceall(line, "%FLAG ", "") + " \n";
          starting = false;
          atomname = false;
          continue;

        } else if (line.find("ATOM_NAME") != string::npos) { //special treatment for atom names activation
          newlines << "END\n" + utils::StringOps::replaceall(line, "%FLAG ", "") + " \n";
          starting = false;
          atomname = true;
          continue;

        } else {
          newlines << "END\n" + utils::StringOps::replaceall(line, "%FLAG ", "") + " \n";
          atomname = false;
          continue;
        }

      } else if (line.find("%") != string::npos) // ignore the other comments
        continue;

      else { //normal line
        if (atomname) { //special case for ATOM_NAME (format: 20a4)
          for (int i = 0; i < line.length() / 4 ; i++) {
            newlines << " " << line.substr(i * 4, 4) << " "; // not optimal check if "   " is correct
          }

          newlines << " " << " \n";

        } else
          newlines << " " << line + " \n";
      }
    }

    amberTopFile.close();
    newlines << "END \n";    //last line
    //output of new lines
    return newlines.str();

  } catch (const gromos::Exception& e) {
    cerr << e.what() << endl;
    exit(1);
  }
}

void AmberTopology::init() {
  if (!stream()) {
    ostringstream msg;
    msg << "Could not open AMBER topology file." << name();
    throw gromos::Exception("AmberTopology", msg.str());
  }

  // set the constants
  d_gff.setFpepsi(gmath::physConst.get_four_pi_eps_i());
  d_gff.setHbar(gmath::physConst.get_hbar());
  d_gff.setSpdl(gmath::physConst.get_speed_of_light());
  d_gff.setBoltz(gmath::physConst.get_boltzmann());
  // First read the whole file into the map
  vector<string> buffer;
  vector<string>::const_iterator it;

  while (!stream().eof()) {
    this->getblock(buffer);

    if (buffer.size()) {
      //std::cout << buffer[0] << std::endl;
      d_blocks[buffer[0]] = buffer;
      buffer.clear();
    }
  }
}

void AmberTopology::parseFile(LinearTopology& lt, double ljscaling, bool atomic_chargegroups) {
  string s, completeblock;
  vector<string> buffer;
  int tmp, n, i;
  double d;
  // Ignored blocks: TITLE, ATOMIC_NUMBER
  // POINTERS
  int natom, natomtype, nbondh, nbond, nangleh, nangle, ndiheh, ndihe, nres;
  int nbondtype, nangletype, ndihetype;
  buffer.clear();
  buffer = d_blocks["POINTERS"];
  gio::concatenate(buffer.begin() + 1, buffer.end() - 1, completeblock);
  _lineStream.clear();
  _lineStream.str(completeblock);
  _lineStream >> natom >> natomtype >> nbondh >> nbond >> nangleh >> nangle >> ndiheh
              >> ndihe >> tmp >> tmp >> tmp >> nres >> tmp >> tmp >> tmp >> nbondtype
              >> nangletype >> ndihetype; // rest is ignored

  if (_lineStream.fail() || _lineStream.eof()) {
    ostringstream msg;
    msg << "In block POINTERS: missing or wrong value";
    throw gromos::Exception("AmberTopology", msg.str());
  }

  //cerr << natom << " " << natomtype << " " << nres << endl;
  // ATOM_NAME
  buffer.clear();
  buffer = d_blocks["ATOM_NAME"];
  completeblock = "";
  gio::concatenate(buffer.begin() + 1, buffer.end() - 1, completeblock);
  _lineStream.clear();
  _lineStream.str(completeblock);

  for (n = 0; n < natom; ++n) {
    _lineStream >> s;

    if (_lineStream.fail() || _lineStream.eof()) {
      ostringstream msg;
      msg << "In block ATOM_NAME: missing or wrong value";
      throw gromos::Exception("AmberTopology", msg.str());
    }

    lt.addAtom(AtomTopology());
    lt.atoms()[n].setName(s);
  }

  //cerr << "ATOM_NAME block read" << endl;
  // CHARGE
  buffer.clear();
  buffer = d_blocks["CHARGE"];
  completeblock = "";
  gio::concatenate(buffer.begin() + 1, buffer.end() - 1, completeblock);
  _lineStream.clear();
  _lineStream.str(completeblock);

  for (n = 0; n < natom; ++n) {
    _lineStream >> d;

    if (_lineStream.fail() || _lineStream.eof()) {
      ostringstream msg;
      msg << "In block CHARGE: missing or wrong value";
      throw gromos::Exception("AmberTopology", msg.str());
    }

    // units: include already 1/(4*PI*eps_0)
    // change to GROMOS units: e
    // conversion factor as used in GROMACS: 18.222615 //TODO Why not amber value? 18.223
    //d /= 18.222615;
    d /= 18.2223; //match parmed for better Gromacs comparison
    d = (int)(10000 * d) / 10000.;
    lt.atoms()[n].setCharge(d);
  }

  //cerr << "CHARGE block read" << endl;
  // MASS
  buffer.clear();
  buffer = d_blocks["MASS"];
  completeblock = "";
  concatenate(buffer.begin() + 1, buffer.end() - 1, completeblock);
  _lineStream.clear();
  _lineStream.str(completeblock);

  for (n = 0; n < natom; ++n) {
    _lineStream >> d;

    if (_lineStream.fail() || _lineStream.eof()) {
      ostringstream msg;
      msg << "In block MASS: missing or wrong value";
      throw gromos::Exception("AmberTopology", msg.str());
    }

    lt.atoms()[n].setMass(d);
  }

  //cerr << "MASS block read" << endl;
  // ATOM_TYPE_INDEX
  vector<int> atom_type_indices(natom);
  buffer.clear();
  buffer = d_blocks["ATOM_TYPE_INDEX"];
  completeblock = "";
  gio::concatenate(buffer.begin() + 1, buffer.end() - 1, completeblock);
  _lineStream.clear();
  _lineStream.str(completeblock);

  for (n = 0; n < natom; ++n) {
    _lineStream >> i;

    if (_lineStream.fail() || _lineStream.eof()) {
      ostringstream msg;
      msg << "In block ATOM_TYPE_INDEX: missing or wrong value";
      throw gromos::Exception("AmberTopology", msg.str());
    }

    lt.atoms()[n].setIac(--i);
    atom_type_indices[n] = i;
  }

  //cerr << "ATOM_TYPE_INDEX block read" << endl;
  // NUMBER_EXCLUDED_ATOMS ignored (not necessary)
  // NONBONDED_PARM_INDEX ignored (not necessary)
  // RESIDUE_LABEL
  buffer.clear();
  buffer = d_blocks["RESIDUE_LABEL"];
  completeblock = "";
  concatenate(buffer.begin() + 1, buffer.end() - 1, completeblock);
  _lineStream.clear();
  _lineStream.str(completeblock);

  for (n = 0; n < nres; ++n) {
    _lineStream >> s;

    if (_lineStream.fail() || _lineStream.eof()) {
      ostringstream msg;
      msg << "In block RESIDUE_LABEL: missing or wrong value";
      throw gromos::Exception("AmberTopology", msg.str());
    }

    lt.setResName(n, s);
  }

  //cerr << "RESIDUE_LABEL block read" << endl;
  // RESIDUE_POINTER
  vector<int> respointer;
  buffer.clear();
  buffer = d_blocks["RESIDUE_POINTER"];
  completeblock = "";
  concatenate(buffer.begin() + 1, buffer.end() - 1, completeblock);
  _lineStream.clear();
  _lineStream.str(completeblock);

  for (n = 0; n < nres; ++n) {
    _lineStream >> i;

    if (_lineStream.fail() || _lineStream.eof()) {
      ostringstream msg;
      msg << "In block RESIDUE_POINTER: missing or wrong value";
      throw gromos::Exception("AmberTopology", msg.str());
    }

    respointer.push_back(i);
  }

  respointer.push_back(natom + 1);

  for (n = 0; n < nres; ++n) {
    int j;

    for (j = respointer[n]; j < respointer[n + 1] - 1; ++j) {
      lt.setResNum(j - 1, n);
      // set charge group
      if(atomic_chargegroups) // each atom belongs to its own chargegroup
        lt.atoms()[j - 1].setChargeGroup(1);
      else // one residue = one charge group
        lt.atoms()[j - 1].setChargeGroup(0);
    }

    // last atom
    lt.setResNum(j - 1, n);
    lt.atoms()[j - 1].setChargeGroup(1);
  }

  //cerr << "RESIDUE_POINTER block read" << endl;
  // BOND_FORCE_CONSTANT
  vector<double> force_constants(nbondtype);
  buffer.clear();
  buffer = d_blocks["BOND_FORCE_CONSTANT"];
  completeblock = "";
  concatenate(buffer.begin() + 1, buffer.end() - 1, completeblock);
  _lineStream.clear();
  _lineStream.str(completeblock);

  for (n = 0; n < nbondtype; ++n) {
    _lineStream >> d;

    if (_lineStream.fail() || _lineStream.eof()) {
      ostringstream msg;
      msg << "In block BOND_FORCE_CONSTANT: missing or wrong value";
      throw gromos::Exception("AmberTopology", msg.str());
    }

    // units: kcal/(mol Angstrom**2)
    // convert into GROMOS units: kJ/(mol nm**2)
    d *= 4.184 * 100;
    force_constants[n] = 2 * d; // in AMBER the factor 0.5 is included
  }

  // BOND_EQUIL_VALUE
  buffer.clear();
  buffer = d_blocks["BOND_EQUIL_VALUE"];
  completeblock = "";
  gio::concatenate(buffer.begin() + 1, buffer.end() - 1, completeblock);
  _lineStream.clear();
  _lineStream.str(completeblock);

  for (n = 0; n < nbondtype; ++n) {
    _lineStream >> d;

    if (_lineStream.fail() || _lineStream.eof()) {
      string msg = "In block BOND_EQUIL_VALUE: missing or wrong value";
      throw gromos::Exception("AmberTopology", msg);
    }

    // units: Angstroem
    // convert into GROMOS units: nm
    d *= 0.1;
    d_gff.addBondType(BondType(n, force_constants[n] / (2 * d * d), force_constants[n], d));
  }

  //cerr << "BONDTYPE block read" << endl;
  // ANGLE_FORCE_CONSTANT
  force_constants.resize(nangletype, 0.0);
  double converter = 180 / gmath::physConst.get_pi();
  buffer.clear();
  buffer = d_blocks["ANGLE_FORCE_CONSTANT"];
  completeblock = "";
  gio::concatenate(buffer.begin() + 1, buffer.end() - 1, completeblock);
  _lineStream.clear();
  _lineStream.str(completeblock);

  for (n = 0; n < nangletype; ++n) {
    _lineStream >> d;

    if (_lineStream.fail() || _lineStream.eof()) {
      std::ostringstream msg;
      msg << "In block ANGLE_FORCE_CONSTANT: missing or wrong value";
      throw gromos::Exception("AmberTopology", msg.str());
    }

    // units: kcal/(mol radian**2)
    // convert into GROMOS units: kJ/(mol degree**2)
    d /= (converter * converter);
    force_constants[n] = d * 4.184 * 2; // in AMBER the factor 0.5 is included
  }

  // ANGLE_EQUIL_VALUE
  buffer.clear();
  buffer = d_blocks["ANGLE_EQUIL_VALUE"];
  completeblock = "";
  gio::concatenate(buffer.begin() + 1, buffer.end() - 1, completeblock);
  _lineStream.clear();
  _lineStream.str(completeblock);

  for (n = 0; n < nangletype; ++n) {
    _lineStream >> d;

    if (_lineStream.fail() || _lineStream.eof()) {
      std::ostringstream msg;
      msg << "In block ANGLE_EQUIL_VALUE: missing or wrong value";
      throw gromos::Exception("AmberTopology", msg.str());
    }

    // units: radian
    // convert into GROMOS units: degree
    double cosd = cos(d);
    d *= converter;
    // quartic force constant (from GROMOS96 book, p. II-20)
    double ekbt = 2.494323; // in kJ/mol
    double t1 = (d + sqrt(ekbt / force_constants[n])) / converter;
    double t2 = (d - sqrt(ekbt / force_constants[n])) / converter;
    double qfc = 2 * ekbt / ((cos(t1) - cosd) * (cos(t1) - cosd) + (cos(t2) - cosd) * (cos(t2) - cosd));
    d_gff.addAngleType(AngleType(n, qfc, force_constants[n], d));
  }

  //cerr << "ANGLETYPE block read" << endl;
# ifdef Debug_DihedType
  cerr << "Dihedral Type start!" << endl;
# endif
  // DIHEDRAL_FORCE_CONSTANT
  force_constants.resize(ndihetype, 0.0);
  buffer.clear();
  buffer = d_blocks["DIHEDRAL_FORCE_CONSTANT"];
  completeblock = "";
  gio::concatenate(buffer.begin() + 1, buffer.end() - 1, completeblock);
  _lineStream.clear();
  _lineStream.str(completeblock);

  for (n = 0; n < ndihetype; ++n) {
    _lineStream >> d;

    if (_lineStream.fail() || _lineStream.eof()) {
      std::ostringstream msg;
      msg << "In block DIHEDRAL_FORCE_CONSTANT: missing or wrong value";
      throw gromos::Exception("AmberTopology", msg.str());
    }

    // units: kcal/mol
    // convert into GROMOS units: kJ/mol
    // AMBER uses a factor 0.5 in the formula, therefore no factor 2
    force_constants[n] = d * 4.184;
  }

  // DIHEDRAL_PERIODICITY
  std::vector<int> periodicity(ndihetype, 0.0);
  buffer.clear();
  buffer = d_blocks["DIHEDRAL_PERIODICITY"];
  completeblock = "";
  gio::concatenate(buffer.begin() + 1, buffer.end() - 1, completeblock);
  _lineStream.clear();
  _lineStream.str(completeblock);

  for (n = 0; n < ndihetype; ++n) {
    _lineStream >> d;

    if (_lineStream.fail() || _lineStream.eof()) {
      std::ostringstream msg;
      msg << "In block DIHEDRAL_PERIODICITY: missing or wrong value";
      throw gromos::Exception("AmberTopology", msg.str());
    }

    if (int(d) < 0) {
      std::ostringstream msg;
      msg << "In block DIHEDRAL_PERIODICITY: negative periodicity values are not supported";
      throw gromos::Exception("AmberTopology", msg.str());
    }

    periodicity[n] = int(d);
  }

  // DIHEDRAL_PHASE
  buffer.clear();
  buffer = d_blocks["DIHEDRAL_PHASE"];
  completeblock = "";
  gio::concatenate(buffer.begin() + 1, buffer.end() - 1, completeblock);
  _lineStream.clear();
  _lineStream.str(completeblock);

  for (n = 0; n < ndihetype; ++n) {
    _lineStream >> d;

    if (_lineStream.fail() || _lineStream.eof()) {
      std::ostringstream msg;
      msg << "In block DIHEDRAL_PHASE: missing or wrong value";
      throw gromos::Exception("AmberTopology", msg.str());
    }

    // units: radian
#ifdef Debug_DihedType
    if (force_constants[n] >= 0) {
      cerr << n << "\t"<< force_constants[n] << "\t" << cos(d) << "\t" << d*converter << "\t" << periodicity[n] << "\n";
    } else {
      cerr << n << "\t"<< -force_constants[n] << "\t" << cos(d) << "\t" << (d*converter)+180.0 << "\t" << periodicity[n] << "\n";
    }
#endif
    // If Force constant is negative, flip sign and dephase by 180 (this is to follow gromos convention). 
    if (force_constants[n] >= 0) {
      d_gff.addDihedralType(DihedralType(n, force_constants[n], cos(d), d*converter, periodicity[n]));
    }else {
      d_gff.addDihedralType(DihedralType(n, -force_constants[n], cos(d), (d*converter)+180.0, periodicity[n]));
    }
  }
#ifdef Debug_DihedType
  cerr << "DIHEDRALTYPE block read" << endl;
#endif
  // IMPROPER DIHEDRALS: we just add the GROMOS improper type
  // FIX: This is a hack
  d_gff.addImproperType(ImproperType(0, 0.051, 0.0)); // planar groups
  // SCEE_SCALE_FACTOR ignored, because set in AMBER block in input file
  // SCNB_SCALE_FACTOR ignored, because set by flag
  // LENNARD_JONES_ACOEF
  std::vector<double> c12;
  buffer.clear();
  buffer = d_blocks["LENNARD_JONES_ACOEF"];
  completeblock = "";
  gio::concatenate(buffer.begin() + 1, buffer.end() - 1, completeblock);
  _lineStream.clear();
  _lineStream.str(completeblock);

  for (n = 0; n < (natomtype * (natomtype + 1) / 2); ++n) {
    _lineStream >> d;

    if (_lineStream.fail() || _lineStream.eof()) {
      std::ostringstream msg;
      msg << "In block LENNARD_JONES_ACOEF: missing or wrong value";
      throw gromos::Exception("AmberTopology", msg.str());
    }

    // units: (kcal/mol) Angstroem**12
    // convert into GROMOS units: (kJ/mol) nm**12
    c12.push_back(d * 4.184 * 1E-12);
  }

  // LENNARD_JONES_BCOEF
  buffer.clear();
  buffer = d_blocks["LENNARD_JONES_BCOEF"];
  completeblock = "";
  gio::concatenate(buffer.begin() + 1, buffer.end() - 1, completeblock);
  _lineStream.clear();
  _lineStream.str(completeblock);
  n = 0;

  for (int j = 0; j < natomtype; ++j) {
    for (int k = 0; k < j + 1; ++k) {
      _lineStream >> d;

      if (_lineStream.fail() || _lineStream.eof()) {
        std::ostringstream msg;
        msg << "In block LENNARD_JONES_BCOEF: missing or wrong value";
        throw gromos::Exception("AmberTopology", msg.str());
      }

      // units: (kcal/mol) Angstroem**6
      // convert into GROMOS units: (kJ/mol) nm**6
      d *= 4.184 * 1E-06;
      d_gff.setLJType(AtomPair(k, j), LJType(c12[n], d, ljscaling * c12[n], ljscaling * d));
      ++n;
    }
  }

  for (int j = 0; j <= natomtype; j++)
    d_gff.setLJType(AtomPair(j, natomtype + 1), LJType(0, 0, 0, 0));

  //cerr << "LJ block read" << endl;
  // adjacency matrix: needed to get the correct order of
  // atoms for the improper dihedrals
  std::map<std::pair<int, int>, int> adjlist;
  // BONDS_INC_HYDROGEN
  buffer.clear();
  buffer = d_blocks["BONDS_INC_HYDROGEN"];
  completeblock = "";
  gio::concatenate(buffer.begin() + 1, buffer.end() - 1, completeblock);
  _lineStream.clear();
  _lineStream.str(completeblock);
  int ib, jb;

  for (n = 0; n < nbondh; ++n) {
    _lineStream >> ib >> jb >> i;

    if (_lineStream.fail() || _lineStream.eof()) {
      std::ostringstream msg;
      msg << "In block BONDS_INC_HYDROGEN: missing or wrong value";
      throw gromos::Exception("AmberTopology", msg.str());
    }

    ib /= 3;
    jb /= 3;

    if (jb < ib) { // swap
      int tmp = ib;
      ib = jb;
      jb = tmp;
    }

    Bond bond(ib, jb);
    bond.setType(--i);
    lt.addBond(bond);
    adjlist[std::make_pair(ib, jb)] = 1;

    // set 1,2-exclusion
    if (lt.atoms()[ib].exclusion().size() == 0) {
      Exclusion* e;
      e = new Exclusion;
      e->insert(jb);
      lt.atoms()[ib].setExclusion(*e);
      delete e;

    } else
      lt.atoms()[ib].exclusion().insert(jb);
  }

  // BONDS_WITHOUT_HYDROGEN
  buffer.clear();
  buffer = d_blocks["BONDS_WITHOUT_HYDROGEN"];
  completeblock = "";
  gio::concatenate(buffer.begin() + 1, buffer.end() - 1, completeblock);
  _lineStream.clear();
  _lineStream.str(completeblock);

  for (n = 0; n < nbond; ++n) {
    _lineStream >> ib >> jb >> i;

    if (_lineStream.fail() || _lineStream.eof()) {
      std::ostringstream msg;
      msg << "In block BONDS_WITHOUT_HYDROGEN: missing or wrong value";
      throw gromos::Exception("AmberTopology", msg.str());
    }

    ib /= 3;
    jb /= 3;

    if (jb < ib) { // swap
      int tmp = ib;
      ib = jb;
      jb = tmp;
    }

    Bond bond(ib, jb);
    bond.setType(--i);
    lt.addBond(bond);
    adjlist[std::make_pair(ib, jb)] = 1;

    // set 1,2-exclusion
    if (lt.atoms()[ib].exclusion().size() == 0) {
      Exclusion* e;
      e = new Exclusion;
      e->insert(jb);
      lt.atoms()[ib].setExclusion(*e);
      delete e;

    } else
      lt.atoms()[ib].exclusion().insert(jb);
  }

  //cerr << "BONDS block read" << endl;
  // ANGLES_INC_HYDROGEN
  buffer.clear();
  buffer = d_blocks["ANGLES_INC_HYDROGEN"];
  completeblock = "";
  gio::concatenate(buffer.begin() + 1, buffer.end() - 1, completeblock);
  _lineStream.clear();
  _lineStream.str(completeblock);
  int it, jt, kt;

  for (n = 0; n < nangleh; ++n) {
    _lineStream >> it >> jt >> kt >> i;

    if (_lineStream.fail() || _lineStream.eof()) {
      std::ostringstream msg;
      msg << "In block ANGLES_INC_HYDROGEN: missing or wrong value";
      throw gromos::Exception("AmberTopology", msg.str());
    }

    it /= 3;
    kt /= 3;

    if (kt < it) {
      int tmp = it;
      it = kt;
      kt = tmp;
    }

    Angle angle(it, jt / 3, kt);
    angle.setType(--i);
    lt.addAngle(angle);

    // set 1,3-exclusion
    if (lt.atoms()[it].exclusion().size() == 0) {
      Exclusion* e;
      e = new Exclusion;
      e->insert(kt);
      lt.atoms()[it].setExclusion(*e);
      delete e;

    } else
      lt.atoms()[it].exclusion().insert(kt);
  }

  // ANGLES_WITHOUT_HYDROGEN
  buffer.clear();
  buffer = d_blocks["ANGLES_WITHOUT_HYDROGEN"];
  completeblock = "";
  gio::concatenate(buffer.begin() + 1, buffer.end() - 1, completeblock);
  _lineStream.clear();
  _lineStream.str(completeblock);

  for (n = 0; n < nangle; ++n) {
    _lineStream >> it >> jt >> kt >> i;

    if (_lineStream.fail() || _lineStream.eof()) {
      std::ostringstream msg;
      msg << "In block ANGLES_WITHOUT_HYDROGEN: missing or wrong value";
      throw gromos::Exception("AmberTopology", msg.str());
    }

    it /= 3;
    kt /= 3;

    if (kt < it) {
      int tmp = it;
      it = kt;
      kt = tmp;
    }

    Angle angle(it, jt / 3, kt);
    angle.setType(--i);
    lt.addAngle(angle);

    // set 1,3-exclusion
    if (lt.atoms()[it].exclusion().size() == 0) {
      Exclusion* e;
      e = new Exclusion;
      e->insert(kt);
      lt.atoms()[it].setExclusion(*e);
      delete e;

    } else
      lt.atoms()[it].exclusion().insert(kt);
  }

  //cerr << "ANGLES block read" << endl;
#ifdef DebugDihedH
  cerr << "Start Debug DihedralH \n";
#endif
  // DIHEDRALS_INC_HYDROGEN
  buffer.clear();
  buffer = d_blocks["DIHEDRALS_INC_HYDROGEN"];
  completeblock = "";
  gio::concatenate(buffer.begin() + 1, buffer.end() - 1, completeblock);
  _lineStream.clear();
  _lineStream.str(completeblock);
  int ip, jp, kp, lp;

  for (n = 0; n < ndiheh; ++n) {
    _lineStream >> ip >> jp >> kp >> lp >> i;
#ifdef DebugDihedH
    cerr << ip << "\t" << jp << "\t" << kp << "\t" << lp << "\t" << i << "\n";
#endif

    if (_lineStream.fail() || _lineStream.eof()) {
      std::ostringstream msg;
      msg << "In block DIHEDRALS_INC_HYDROGEN: missing or wrong value";
      throw gromos::Exception("AmberTopology", msg.str());
    }

    ip /= 3;
    jp /= 3;
    kp /= 3;
    lp /= 3;
#ifdef DebugDihedH
    cerr << ip << "\t" << jp << "\t" << kp << "\t" << lp << "\t" << i << "\n";
#endif
    //lp = abs(lp);

    if (lp < 0) { // improper
      // central atom is in AMBER at the second or third position
      // in GROMOS it's at the first position
      if (adjlist.count(std::make_pair(ip, jp)) > 0 ||
          adjlist.count(std::make_pair(jp, ip)) > 0) { // second position
        Dihedral dihedral(abs(ip), abs(jp), abs(kp),
                          abs(lp));    //Todo: AmberImpropers as Dihedals - bschroed //added abs GK
        dihedral.setType(i - 1);                        //Todo: AmberImpropers as Dihedals - bschroed
        lt.addDihedral(dihedral);                       //Todo: AmberImpropers as Dihedals - bschroed
        //Improper improper(jp, ip, abs(kp), abs(lp));
        //improper.setType(0); // FIX: this is a hack!!!
        //lt.addImproper(improper);

      } else { // third position
        // Use impropers as Dihedrals (amberstyle)
        Dihedral dihedral(abs(ip), abs(jp), abs(kp),
                          abs(lp));    //Todo: AmberImpropers as Dihedals - bschroed //added abs GK
        dihedral.setType(i - 1);                        //Todo: AmberImpropers as Dihedals - bschroed
        lt.addDihedral(dihedral);                       //Todo: AmberImpropers as Dihedals - bschroed
        //use impropers but yet only with one type (gromosstyle)
        //Improper improper(abs(kp), ip, jp, abs(lp));
        //improper.setType(0); // FIX: this is a hack!!!
        //lt.addImproper(improper);
      }

    } else { // proper
      //if (lp < ip) { // i should always be the smaller
      if (abs(kp) < abs(
            jp)) { // GK: According to the manual j should always be the smaller than k (however, then no dihe are written )
        int tmp = ip;
        ip = lp;
        lp = tmp;
        tmp = jp;
        jp = kp;
        kp = tmp;
      }

      if (force_constants[i - 1] > 0) { // we only add torsions with a non-zero force constant
        Dihedral dihedral(abs(ip), abs(jp), abs(kp), abs(lp)); //added abs GK
        dihedral.setType(i - 1);
        lt.addDihedral(dihedral);
      }

      if (!lt.atoms()[min(abs(ip), abs(lp))].exclusion().contains(max(abs(ip),
                                                                      abs(lp))) // using min/max of ip and lp here since the dihedrals have to be ordered based on j and k !  GK
          && !lt.atoms()[min(abs(ip), abs(lp))].exclusion14().contains(max(abs(ip),
                                                                           abs(lp)))) { // set 1,4-exclusion
        if (lt.atoms()[min(abs(ip), abs(lp))].exclusion14().size() == 0) {
          Exclusion* e;
          e = new Exclusion;
          e->insert(max(abs(ip), abs(lp)));
          lt.atoms()[min(abs(ip), abs(lp))].setExclusion14(*e);
          delete e;

        } else
          lt.atoms()[min(abs(ip), abs(lp))].exclusion14().insert(max(abs(ip), abs(lp)));
      }
    }
  }

#ifdef DebugDihedH
  cerr << "END Debug DihedralH \n";
#endif
#ifdef DebugDihedwoH
  cerr << "Start Debug DihedralwoH \n";
#endif
  // DIHEDRALS_WITHOUT_HYDROGEN
  buffer.clear();
  buffer = d_blocks["DIHEDRALS_WITHOUT_HYDROGEN"];
  completeblock = "";
  gio::concatenate(buffer.begin() + 1, buffer.end() - 1, completeblock);
  _lineStream.clear();
  _lineStream.str(completeblock);

  for (n = 0; n < ndihe; ++n) {
    _lineStream >> ip >> jp >> kp >> lp >> i;
#ifdef DebugDihedwoH
    cerr << ip << "\t" << jp << "\t" << kp << "\t" << lp << "\t" << i << "\n";
#endif

    if (_lineStream.fail() || _lineStream.eof()) {
      std::ostringstream msg;
      msg << "In block DIHEDRALS_WITHOUT_HYDROGEN: missing or wrong value";
      throw gromos::Exception("AmberTopology", msg.str());
    }

    ip = (ip / 3);
    jp = (jp / 3);
    kp = (kp / 3);
    lp = (lp / 3);

    //lp = abs(lp);
    if (lp < 0) { // improper
      // central atom is in AMBER at the second or third position
      // in GROMOS it's at the first position
      if (adjlist.count(std::make_pair(ip, jp)) > 0 ||
          adjlist.count(std::make_pair(jp, ip)) > 0) { // second position
        Dihedral dihedral(abs(ip), abs(jp), abs(kp),
                          abs(lp));  //Todo: AmberImpropers as Dihedals - bschroed //added abs GK
        dihedral.setType(i - 1);  //Todo: AmberImpropers as Dihedals - bschroed
        lt.addDihedral(dihedral); //Todo: AmberImpropers as Dihedals - bschroed
        //Improper improper(jp, ip, abs(kp), abs(lp));
        //improper.setType(0); // FIX: this is a hack!!!
        //lt.addImproper(improper);

      } else { // third position
        Dihedral dihedral(abs(ip), abs(jp), abs(kp),
                          abs(lp));  //Todo: AmberImpropers as Dihedals - bschroed //added abs GK
        dihedral.setType(i - 1);  //Todo: AmberImpropers as Dihedals - bschroed
        lt.addDihedral(dihedral); //Todo: AmberImpropers as Dihedals - bschroed
        //Improper improper(abs(kp), ip, jp, abs(lp));
        //improper.setType(0); // FIX: this is a hack!!!
        //lt.addImproper(improper);
      }

    } else { // proper
      //if (lp < ip) { // i should always be the smaller
      if (abs(kp) < abs(
            jp)) { // GK: According to the manual j should always be the smaller than k (however, then no dihe are written )
        int tmp = ip;
        ip = lp;
        lp = tmp;
        tmp = jp;
        jp = kp;
        kp = tmp;
      }

      if (force_constants[i - 1] > 0) { // we only add torsions with a non-zero force constant
        Dihedral dihedral(abs(ip), abs(jp), abs(kp), abs(lp)); //added abs GK
        dihedral.setType(i - 1);
        lt.addDihedral(dihedral);
      }

      if (!lt.atoms()[min(abs(ip), abs(lp))].exclusion().contains(max(abs(ip),
                                                                      abs(lp))) // using min/max of ip and lp here since the dihedrals have to be ordered based on j and k !  GK
          && !lt.atoms()[min(abs(ip), abs(lp))].exclusion14().contains(max(abs(ip),
                                                                           abs(lp)))) { // set 1,4-exclusion
        if (lt.atoms()[min(abs(ip), abs(lp))].exclusion14().size() == 0) {
          Exclusion* e;
          e = new Exclusion;
          e->insert(max(abs(ip), abs(lp)));
          lt.atoms()[min(abs(ip), abs(lp))].setExclusion14(*e);
          delete e;

        } else
          lt.atoms()[min(abs(ip), abs(lp))].exclusion14().insert(max(abs(ip), abs(lp)));
      }
    }
  }

  //cerr << "DIHEDRALS block read" << endl;
#ifdef DebugDihedwoH
  cerr << "END Debug DihedralwoH \n";
#endif
  // AMBER_ATOM_TYPE
  std::vector<std::vector<string>> amber_atom_types(natomtype + 1);
  buffer.clear();
  buffer = d_blocks["AMBER_ATOM_TYPE"];
  completeblock = "";
  gio::concatenate(buffer.begin() + 1, buffer.end() - 1, completeblock);
  _lineStream.clear();
  _lineStream.str(completeblock);

  // GK: Initialize atom type array with some dummy atom types to avoid segfaults
  for (i = 0; i < natomtype; ++i)
    amber_atom_types[i].push_back("XYZA");

  for (n = 0; n < natom; ++n) {
    _lineStream >> s;

    if (_lineStream.fail() || _lineStream.eof()) {
      std::ostringstream msg;
      msg << "In block AMBER_ATOM_TYPE: missing or wrong value";
      throw gromos::Exception("AmberTopology", msg.str());
    }

    bool ovwr = true ;  // GK Permission to overwrite dummy type

    for (int i = 0; i < n ;
         i++) { // GK Go through list of previous atom types to avoid using the same name for two types
      if (amber_atom_types[atom_type_indices[i]][0] == s) {
        ovwr = false;
        break;
      };
    } //

    if (ovwr)
      amber_atom_types[atom_type_indices[n]][0] = s;; // GK
  }

  int xcounter = 0 ; // GK count how many dummy types are required

  for (n = 0; n < natomtype; ++n) {
    if (amber_atom_types[n][0] == "XYZA") {
      amber_atom_types[n][0] = "X" + std::to_string(xcounter);
      xcounter++;
    }

    d_gff.addAtomTypeName(amber_atom_types[n][0]); // we just take the first one
  }

  d_gff.addAtomTypeName("DUM");
  natomtype++;
}
}
