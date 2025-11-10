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

// gio_InBuildingBlock.cc

#include "InLinkingBlock.h"

#include <cassert>
#include <cmath>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

#include "Ginstream.h"
#include "../gcore/BbLink.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/Bond.h"
#include "../gcore/Angle.h"
#include "../gcore/Dihedral.h"
#include "../gcore/Improper.h"
#include "../gcore/Exclusion.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/BbSolute.h"
#include "../gcore/AtomTopology.h"

using namespace gcore;
using gio::InLinkingBlock;

// Implementation class

class gio::InLinkingBlock_i : public gio::Ginstream {
  friend class gio::InLinkingBlock;
  gcore::BbLink d_lnk;
  std::vector<gcore::BbLink> d_lnks;
  std::string d_ff;
  /**
   * The init function reads in the whole file and stores
   * the information in the BuildingBlock
   */
  void init();
  /**
   * A function to read in MTBUILDBLSOLUTE blocks
   */
  void readLink(std::vector<std::string> & buffer);
  void readForceField(std::vector<std::string> &buffer);

  void readAtoms(BbLink &bb, std::string bname, std::string resname);
  void readBonds(BbLink &bb, std::string bname, std::string resname);
  void readAngles(BbLink &bb, std::string bname, std::string resname);
  void readImpropers(BbLink &bb, std::string bname, std::string resname);
  void readDihedrals(BbLink &bb, std::string bname, std::string resname);

  InLinkingBlock_i(std::string &s) : d_lnks() {
    this->open(s);
    this->init();
  }
};

// Constructors

InLinkingBlock::InLinkingBlock(std::string name) {
  d_this = new InLinkingBlock_i(name);
}

InLinkingBlock::~InLinkingBlock() {
  delete d_this;
}

const std::string InLinkingBlock::title()const {
  return d_this->title();
}

const std::vector<gcore::BbLink> &InLinkingBlock::links()const {
  return d_this->d_lnks;
}

void gio::InLinkingBlock_i::init() {
  if (!stream())
    throw InLinkingBlock::Exception("Could not open Linking block file "
          + name());

  std::vector<std::string> buffer;

  while (!stream().eof()) {
    getblock(buffer);
    if (buffer.size()) {
      if (buffer[0] == "MTBUILDBLLINK") readLink(buffer);
      else if (buffer[0] == "FORCEFIELD") readForceField(buffer);
      else
	throw InLinkingBlock::Exception("Don't know how to handle block"
	+buffer[0]);
    }
  }
}


void gio::InLinkingBlock_i::readForceField(std::vector<std::string> &buffer) {
  // This block contains one line with a force field code
  if (buffer.size() != 3)
    throw InLinkingBlock::Exception("LinkingBlock file " + name() +
          " is corrupted. FORCEFIELD block should have only one line");
  if (buffer[buffer.size() - 1].find("END") != 0)
    throw InLinkingBlock::Exception("LinkingBlock file " + name() +
          " is corrupted. No END in FORCEFIELD"
          " block. Got\n"
          + buffer[buffer.size() - 1]);

  d_ff = buffer[1];
}

void gio::InLinkingBlock_i::readLink(std::vector<std::string> &buffer) {
  if (buffer.size() < 3)
    throw InLinkingBlock::Exception("Linking block file " + name()
          + " is corrupted. Empty MTBUILDBLLINK block.");
  if (buffer[buffer.size() - 1].find("END") != 0)
    throw InLinkingBlock::Exception("LinkingBlock file " + name() +
          " is corrupted. No END in"
          " MTBUILDBLLINK block (" + buffer[1] +
          "). Got\n" + buffer[buffer.size() - 1]);

  // generic variables
  std::string resname, s;

  std::string block;
  BbLink bb;
  gio::concatenate(buffer.begin() + 1, buffer.end() - 1, block);

  _lineStream.clear();
  _lineStream.str(block);
  _lineStream >> resname;

  // std::cerr << "==================================================" << std::endl;
  // std::cerr << "BLOCK " << resname << std::endl;
  // std::cerr << "==================================================" << std::endl;
  if (_lineStream.fail())
    throw InLinkingBlock::Exception("Bad line in MTBUILDBLLINK block:\n"
          + block + "Trying to read a linker name");
  bb.setResName(resname);

  // Atoms
  readAtoms(bb, block, resname);

  // Bonds
  readBonds(bb, block, resname);

  // Angles
  readAngles(bb, block, resname);

  // Impropers
  readImpropers(bb, block, resname);

  // Dihedrals
  readDihedrals(bb, block, resname);

  _lineStream >> s;
  if (!_lineStream.eof())
    throw InLinkingBlock::Exception("Bad line in MTBUILDBLLINK block "
          + resname + ".\nTrailing data after dihedrals:" + s);

  d_lnks.push_back(bb);
}

void gio::InLinkingBlock_i::readAtoms(BbLink &bb, std::string bname, std::string resname) {
  double d[3];
  int i[5], na;
  std::string s;

  _lineStream >> na;
  if (_lineStream.fail())
    throw InLinkingBlock::Exception("Bad line in " + bname + " block "
          + resname + ":\nTrying to read NMAT");

  for (int j = 0; j < na; j++) {
    AtomTopology at;
    _lineStream >> i[0];
    if (i[0] != j + 1) {
      std::ostringstream os;
      os << "Atom numbers in " + bname + " block "
              << resname
              << " are not sequential\n"
              << " got " << i[0] << " - expected " << j + 1;

      throw InLinkingBlock::Exception(os.str());
    }
    // read the residue identifier
    _lineStream >> i[0];
    bb.setLinkRes(j,i[0]-1);
    _lineStream >> s;
    at.setName(s);
    _lineStream >> i[0] >> d[0] >> d[1] >> i[1];
    if (_lineStream.fail()) {
      std::ostringstream os;
      os << "Bad line in " + bname + " block " << resname
              << ".\nTrying to read atom number " << j + 1;
      throw InLinkingBlock::Exception(os.str());
    }
    at.setIac(--i[0]);
    // WARNING: in the building block we use the mass code, 
    //          in the AtomTopology we usually have real masses
    at.setMass(d[0]);
    at.setCharge(d[1]);
    at.setChargeGroup(i[1]);

    // The trailing atoms do not have exclusions specified.
    // these are specified as the preceding exclusions of the 
    // next residue
    Exclusion e;
    _lineStream >> i[0];
    for (int k = 0; k < i[0]; k++) {
      _lineStream >> i[1];
      e.insert(--i[1]);
    }
    if (_lineStream.fail()) {
      std::ostringstream os;
      os << "Bad line in " + bname + " block " << resname
              << ".\nTrying to read exclusions of atom " << j + 1;
      throw InLinkingBlock::Exception(os.str());
    }
    at.setExclusion(e);
    bb.addAtom(at);
  }
}

void gio::InLinkingBlock_i::readBonds(BbLink &bb, std::string bname, std::string resname) {
  int i[3], num;
  _lineStream >> num;

  if (_lineStream.fail())
    throw InLinkingBlock::Exception("Bad line in " + bname + " block "
          + resname + ".\nTrying to read number of bonds.");
  for (int j = 0; j < num; j++) {
    _lineStream >> i[0] >> i[1] >> i[2];
    if (_lineStream.fail()) {
      std::ostringstream os;
      os << "Bad line in " + bname + " block " << resname
              << ".\nTrying to read " << num << " bonds";
      throw InLinkingBlock::Exception(os.str());
    }
    Bond bond(--i[0], --i[1]);
    bond.setType(--i[2]);
    bb.addBond(bond);
  }
}

void gio::InLinkingBlock_i::readAngles(BbLink &bb, std::string bname, std::string resname) {
  int i[4], num;
  _lineStream >> num;

  if (_lineStream.fail())
    throw InLinkingBlock::Exception("Bad line in " + bname + " block "
          + resname + ".\nTrying to read number of angles.");
  for (int j = 0; j < num; j++) {
    _lineStream >> i[0] >> i[1] >> i[2] >> i[3];
    if (_lineStream.fail()) {
      std::ostringstream os;
      os << "Bad line in " + bname + " block " << resname
              << ".\nTrying to read " << num << " angles";
      throw InLinkingBlock::Exception(os.str());
    }
    Angle angle(--i[0], --i[1], --i[2]);
    angle.setType(--i[3]);
    bb.addAngle(angle);
  }
}

void gio::InLinkingBlock_i::readImpropers(BbLink &bb, std::string bname, std::string resname) {
  int i[5], num;
  _lineStream >> num;

  if (_lineStream.fail())
    throw InLinkingBlock::Exception("Bad line in " + bname + " block "
          + resname + ".\nTrying to read number of impropers.");
  for (int j = 0; j < num; j++) {
    _lineStream >> i[0] >> i[1] >> i[2] >> i[3] >> i[4];
    if (_lineStream.fail()) {
      std::ostringstream os;
      os << "Bad line in " + bname + " block " << resname
              << ".\nTrying to read " << num << " impropers";
      throw InLinkingBlock::Exception(os.str());
    }
    Improper improper(--i[0], --i[1], --i[2], --i[3]);
    improper.setType(--i[4]);
    bb.addImproper(improper);
  }
}

void gio::InLinkingBlock_i::readDihedrals(BbLink &bb, std::string bname, std::string resname) {
  int i[5], num;
  _lineStream >> num;

  if (_lineStream.fail())
    throw InLinkingBlock::Exception("Bad line in " + bname + " block "
          + resname + ".\nTrying to read number of dihedrals.");
  for (int j = 0; j < num; j++) {
    _lineStream >> i[0] >> i[1] >> i[2] >> i[3] >> i[4];
    if (_lineStream.fail()) {
      std::ostringstream os;
      os << "Bad line in " + bname + " block " << resname
              << ".\nTrying to read " << num << " dihedrals\n"
              << i[0] << " " << i[1] << " " << i[2] << " " << i[3]
              << " of type " << i[4];
      throw InLinkingBlock::Exception(os.str());
    }
    Dihedral dihedral(--i[0], --i[1], --i[2], --i[3]);
    dihedral.setType(--i[4]);
    bb.addDihedral(dihedral);
  }
}

