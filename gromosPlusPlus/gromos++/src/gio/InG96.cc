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

// gio_InG96.cc
#include "InG96.h"

#include <cassert>
#include <cstddef>
#include <iostream>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "Ginstream.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gmath/Vec.h"
#include "../gmath/Matrix.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/Box.h"
#include "../gcore/Remd.h"
#include "../utils/groTime.h"

enum blocktype {
  titleblock, timestep,
  positionred, position,
  velocityred, velocity, cosdisplacements,
  box, triclinicbox, genbox, gmxbox, remd, latticeshifts
};

typedef std::map<std::string, blocktype>::value_type BT;
// Define class variable BT (block_types)
const BT blocktypes[] = {BT("TITLE", titleblock),
  BT("TIMESTEP", timestep),
  BT("POSITIONRED", positionred),
  BT("POSITION", position),
  BT("LATTICESHIFTS", latticeshifts),
  BT("VELOCITYRED", velocityred),
  BT("VELOCITY", velocity),
  BT("COSDISPLACEMENTS", cosdisplacements),
  BT("BOX", box),
  BT("TRICLINICBOX", triclinicbox),
  BT("GENBOX", genbox),
  BT("GMXBOX", gmxbox),
  BT("REMD", remd)};

static std::map<std::string, blocktype> BLOCKTYPE(blocktypes, blocktypes + 13);

using gio::InG96;
using namespace gcore;

class gio::InG96_i : public gio::Ginstream {
  friend class gio::InG96;

  std::string d_current;
  int d_switch;

  int d_skip;
  int d_stride;

  // the TIMESTEP block information
  bool d_time_read;
  int d_step;
  double d_time;

  InG96_i(const std::string &name, int skip, int stride)
  : d_current(),
  d_switch(),
  d_skip(skip),
  d_stride(stride),
  d_time_read(false),
  d_step(0),
  d_time(0.0) {
    open(name);
    getline(d_current);
    d_switch = 0;
  }

  ~InG96_i() {
  }

  // method
  void readTimestep();
  void readPosition(gcore::System &sys);
  void readLatticeshifts(gcore::System &sys);
  void readVelocity(gcore::System &sys);
  void initializeBfactors(gcore::System &sys);
  void readCosDisplacements(gcore::System &sys);
  void readTriclinicbox(gcore::System &sys);
  void readGmxbox(gcore::System &sys);
  void readBox(gcore::System &sys);
  void readGenbox(gcore::System &sys);
  void readRemd(gcore::System &sys);
};

// Constructors

InG96::InG96(const std::string &name, int skip, int stride)
: d_skip(skip),
d_stride(stride),
d_stride_eof(false) {
  d_this = new InG96_i(name, skip, stride);
}

InG96::InG96(int skip, int stride)
: d_this(NULL),
d_skip(skip),
d_stride(stride),
d_stride_eof(false) {
}

InG96::~InG96() {
  if (d_this)delete d_this;
}

void InG96::open(const std::string &name) {
  if (d_this) {
    // recover skip and stride
    d_skip = d_this->d_skip;
    d_stride = d_this->d_stride;

    delete d_this;
  }

  d_this = new InG96_i(name, d_skip, d_stride);
}

void InG96::close() {
  if (d_this) {
    d_this->close();

    d_skip = d_this->d_skip;
    d_stride = d_this->d_stride;

    delete d_this;
    d_this = NULL;
  }
}

void InG96::select(const std::string &thing) {
  if (!d_this) {
    throw InG96::Exception("select must be called after open!");
  }

  if (thing == "ALL") {
    d_this->d_switch = 1;
  } else if (thing == "SOLVENT") {
    d_this->d_switch = 2;
  }
  else if (thing == "SOLUTE") {
    d_this->d_switch = 0;
  } else {
    d_this->d_switch = 0;
  }
}

bool InG96::eof()const {
  if (!d_this) {
    throw InG96::Exception("eof, but no open file");
  }

  return d_this->stream().eof();
}

std::string InG96::title()const {
  if (!d_this) {
    throw InG96::Exception("no open file");
  }
  return d_this->title();
}

void gio::InG96_i::readTimestep() {
  std::vector<std::string> buffer;
  getblock(buffer);
  if (buffer[buffer.size() - 1].find("END") != 0) {
    throw InG96::Exception("Coordinate file " + name() +
            " is corrupted. No END in TIMESTEP"
            " block. Got\n"
            + buffer[buffer.size() - 1]);
  }
  d_time_read = false;
  std::vector<std::string>::iterator it = buffer.begin();
  std::istringstream line(*it);

  line >> d_step >> d_time;
  if (line.fail()) {
    throw InG96::Exception("Coordinate file " + name() +
            " is corrupted. Bad line in TIMESTEP block. Got\n" +
            *it);
  }
  // we read the system
  d_time_read = true;
}

void gio::InG96_i::readPosition(gcore::System &sys) {
  std::vector<std::string> buffer;
  std::vector<std::string>::iterator it;
  getblock(buffer);
  if (buffer[buffer.size() - 1].find("END") != 0)
    throw InG96::Exception("Coordinate file " + name() +
          " is corrupted. No END in POSITION"
          " block. Got\n"
          + buffer[buffer.size() - 1]);

  it = buffer.begin();
  std::string dummy;
  int begin = 0;
  if (d_current == "POSITION") begin = 24;
  int na = 0;
  for (int m = 0; m < sys.numMolecules(); m++) {
    if (!sys.mol(m).numPos()) {
      sys.mol(m).initPos();
    }
    na += sys.mol(m).numAtoms();
  }

  // solute?
  if (d_switch < 2) {
    for (int m = 0; m < sys.numMolecules(); m++) {
      for (int a = 0; a < sys.mol(m).numAtoms(); a++, ++it) {
        if (unsigned (begin) >= (*it).size()) {
          std::ostringstream os;
          os << "Coordinate file " << name() << " corrupted.\n"
                  << "Failed to read coordinates from line \n"
                  << *it
                  << "\nfrom POSITION or POSITIONRED block";
          throw InG96::Exception(os.str());
        }

        _lineStream.clear();
        _lineStream.str((*it).substr(begin, (*it).size()));
        _lineStream >> sys.mol(m).pos(a)[0]
                >> sys.mol(m).pos(a)[1]
                >> sys.mol(m).pos(a)[2];


        if (_lineStream.fail()) {
          std::ostringstream os;
          os << "Coordinate file " << name() << " corrupted.\n"
                  << "Failed to read " << na << " solute coordinates"
                  << " from POSITION or POSITIONRED block";
          throw InG96::Exception(os.str());
        }
      }

    }
  } else it += na;

  // Solvent?
  if (d_switch > 0) {
    sys.sol(0).setNumPos(0);
    gmath::Vec v;

    for (; it != buffer.end() - 1; ++it) {
      if (unsigned(begin) >= (*it).size()) {
        std::ostringstream os;
        os << "Coordinate file " << name() << " corrupted.\n"
                << "Failed to read coordinates from line \n"
                << *it
                << "\nfrom POSITION or POSITIONRED block";
        throw InG96::Exception(os.str());
      }
      _lineStream.clear();
      _lineStream.str((*it).substr(begin, (*it).size()));
      _lineStream >> v[0] >> v[1] >> v[2];

      if (_lineStream.fail()) {
        std::ostringstream os;
        os << "Coordinate file " << name() << " corrupted.\n"
                << "Failed while reading solvent coordinates"
                << " from POSITION or POSITIONRED block";
        throw InG96::Exception(os.str());
      }
      sys.sol(0).addPos(v);
    }

    if (sys.sol(0).numPos() % sys.sol(0).topology().numAtoms() != 0) {
      std::ostringstream os;
      os << "Coordinate file " << name() << " corrupted.\n"
              << "Atom count mismatch while reading solvent coordinates.\n"
              << "Read " << sys.sol(0).numPos() << " coordinates for solvent "
              << "with " << sys.sol(0).topology().numAtoms() << " atoms per molecule\n";
      throw InG96::Exception(os.str());
    }
  }
}

void gio::InG96_i::readLatticeshifts(gcore::System &sys) {
  // mk_script just has to know that there is a block LATTICESHIFTS but
  // does not need the values, so nothing is done here
  // (unless you need the shifts, then implement this function ;-))
}

void gio::InG96_i::initializeBfactors(gcore::System &sys) {
  for (int m = 0; m < sys.numMolecules(); m++) {
    if (!sys.mol(m).numBfac()) {
      sys.mol(m).initBfac();
    }
  }
}

void gio::InG96_i::readVelocity(gcore::System &sys) {
  std::vector<std::string> buffer;
  std::vector<std::string>::iterator it;
  getblock(buffer);
  if (buffer[buffer.size() - 1].find("END") != 0)
    throw InG96::Exception("Coordinate file " + name() +
          " is corrupted. No END in VELOCITY"
          " block. Got\n"
          + buffer[buffer.size() - 1]);

  it = buffer.begin();
  std::string dummy;
  int begin = 0;
  if (d_current == "VELOCITY") begin = 24;
  int na = 0;
  for (int m = 0; m < sys.numMolecules(); m++) {
    if (!sys.mol(m).numVel()) {
      sys.mol(m).initVel();
    }
    na += sys.mol(m).numVel();
  }
  // solute?
  if (d_switch < 2) {
    for (int m = 0; m < sys.numMolecules(); m++) {
      for (int a = 0; a < sys.mol(m).numVel(); a++, ++it) {
        if (unsigned (begin) >= (*it).size()) {
          std::ostringstream os;
          os << "Coordinate file " << name() << " corrupted.\n"
                  << "Failed to read velocities from line \n"
                  << *it
                  << "\nfrom VELOCITY or VELOCITYRED block";
          throw InG96::Exception(os.str());
        }
        _lineStream.clear();
        _lineStream.str((*it).substr(begin, (*it).size()));
        _lineStream >> sys.mol(m).vel(a)[0]
                >> sys.mol(m).vel(a)[1]
                >> sys.mol(m).vel(a)[2];

        if (_lineStream.fail()) {
          std::ostringstream os;
          os << "Coordinate file " << name() << " corrupted.\n"
                  << "Failed to read " << na << " solute velocity coordinates"
                  << " from VELOCITY or VELOCITYRED block";
          throw InG96::Exception(os.str());
        }
      }
    }
  } else it += na;

  // Solvent?
  if (d_switch > 0) {
    sys.sol(0).setNumVel(0);
    gmath::Vec v;

    for (; it != buffer.end() - 1; ++it) {
      if (unsigned (begin) >= (*it).size()) {
        std::ostringstream os;
        os << "Coordinate file " << name() << " corrupted.\n"
                << "Failed to read velocities from line \n"
                << *it
                << "\nfrom VELOCITTY or VELOCITYRED block";
        throw InG96::Exception(os.str());
      }
      _lineStream.clear();
      _lineStream.str((*it).substr(begin, (*it).size()));
      _lineStream >> v[0] >> v[1] >> v[2];
      sys.sol(0).addVel(v);

      if (_lineStream.fail()) {
        std::ostringstream os;
        os << "Coordinate file " << name() << " corrupted.\n"
                << "Failed while reading solvent velocity coordinates"
                << " from VELOCITY or VELOCITYRED block";
        throw InG96::Exception(os.str());
      }
    }
    if (sys.sol(0).numVel() % sys.sol(0).topology().numAtoms() != 0) {
      std::ostringstream os;
      os << "Coordinate file " << name() << " corrupted.\n"
              << "Atom count mismatch while reading solvent velocities.\n"
              << "Read " << sys.sol(0).numVel() << " coordinates for solvent "
              << "with " << sys.sol(0).topology().numAtoms() << " atoms per molecule\n";
      throw InG96::Exception(os.str());
    }
  }
}

void gio::InG96_i::readCosDisplacements(gcore::System &sys) {
  std::vector<std::string> buffer;
  std::vector<std::string>::iterator it;
  getblock(buffer);
  if (buffer[buffer.size() - 1].find("END") != 0)
    throw InG96::Exception("Coordinate file " + name() +
          " is corrupted. No END in COSDISPLACEMENTS"
          " block. Got\n"
          + buffer[buffer.size() - 1]);

  it = buffer.begin();
  int na = 0;
  for (int m = 0; m < sys.numMolecules(); m++) {
    if (!sys.mol(m).numCosDisplacements()) {
      sys.mol(m).initCosDisplacements();
    }
    na += sys.mol(m).numCosDisplacements();
  }
  // solute?
  if (d_switch < 2) {
    for (int m = 0; m < sys.numMolecules(); m++) {
      for (int a = 0; a < sys.mol(m).numCosDisplacements(); a++, ++it) {
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> sys.mol(m).cosDisplacement(a)[0]
                >> sys.mol(m).cosDisplacement(a)[1]
                >> sys.mol(m).cosDisplacement(a)[2];

        if (_lineStream.fail()) {
          std::ostringstream os;
          os << "Coordinate file " << name() << " corrupted.\n"
                  << "Failed to read " << na << " solute COS displacements"
                  << " from COSDISPLACEMENTS block";
          throw InG96::Exception(os.str());
        }
      }
    }
  } else it += na;

  // Solvent?
  if (d_switch > 0) {
    sys.sol(0).setNumCosDisplacements(0);
    gmath::Vec v;

    for (; it != buffer.end() - 1; ++it) {
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> v[0] >> v[1] >> v[2];
      sys.sol(0).addCosDisplacement(v);

      if (_lineStream.fail()) {
        std::ostringstream os;
        os << "Coordinate file " << name() << " corrupted.\n"
                << "Failed while reading solvent COS displacements"
                << " from COSDISPLACEMENTS block.";
        throw InG96::Exception(os.str());
      }
    }
    if (sys.sol(0).numCosDisplacements() % sys.sol(0).topology().numAtoms() != 0) {
      std::ostringstream os;
      os << "Coordinate file " << name() << " corrupted.\n"
              << "Atom count mismatch while reading solvent COS displacements.\n"
              << "Read " << sys.sol(0).numCosDisplacements() << " coordinates for solvent "
              << "with " << sys.sol(0).topology().numAtoms() << " atoms per molecule\n";
      throw InG96::Exception(os.str());
    }
  }
}

void gio::InG96_i::readBox(gcore::System &sys) {
  // std::cerr << "readbox" << std::endl;
  std::vector<std::string> buffer;
  std::vector<std::string>::iterator it;
  getblock(buffer);
  if (buffer[buffer.size() - 1].find("END") != 0)
    throw InG96::Exception("Coordinate file " + name() +
          " is corrupted. No END in BOX"
          " block. Got\n"
          + buffer[buffer.size() - 1]);

  it = buffer.begin();

  _lineStream.clear();
  _lineStream.str(*it);
  _lineStream >> sys.box().K()[0] >> sys.box().L()[1] >> sys.box().M()[2];
  if (_lineStream.fail())
    throw InG96::Exception("Bad line in BOX block:\n" + *it +
          "\nTrying to read three doubles");
}

void gio::InG96_i::readGmxbox(System &sys) {
  std::vector<std::string> buffer;
  std::vector<std::string>::iterator it;

  getblock(buffer);
  if (buffer[buffer.size() - 1].find("END") != 0)
    throw InG96::Exception("Coordinate file " + name() +
          " is corrupted. No END in GMXBOX"
          " block. Got\n"
          + buffer[buffer.size() - 1]);

  std::string s;
  gio::concatenate(buffer.begin(), buffer.end() - 1, s);

  _lineStream.clear();
  _lineStream.str(s);
  sys.box().setNtb(gcore::Box::boxshape_enum(gcore::Box::triclinic));

  gmath::Vec k, l, m;
  _lineStream >> k[0] >> l[1] >> m[2]
          >> k[1] >> k[2] >> l[0]
          >> l[2] >> m[0] >> m[1];
  sys.box().K() = k;
  sys.box().L() = l;
  sys.box().M() = m;

  sys.box().update_triclinic();

  if (_lineStream.fail())
    throw InG96::Exception("Bad line in GMXBOX block:\n" + *it +
          "\nTrying to read nine doubles");

}

void gio::InG96_i::readTriclinicbox(System &sys) {
  std::vector<std::string> buffer;
  std::vector<std::string>::iterator it;

  getblock(buffer);
  if (buffer[buffer.size() - 1].find("END") != 0)
    throw InG96::Exception("Coordinate file " + name() +
          " is corrupted. No END in TRICLINICBOX"
          " block. Got\n"
          + buffer[buffer.size() - 1]);

  std::string s;
  gio::concatenate(buffer.begin(), buffer.end() - 1, s);

  _lineStream.clear();
  _lineStream.str(s);
  int ntb;
  _lineStream >> ntb;

  // GromosXX truncated octahedral box:  3
  // Gromos++ truncated octahedral box: -1
  // yes, i know!
  //if (ntb == -1) ntb = 3;

  sys.box().setNtb(gcore::Box::boxshape_enum(ntb));

  gmath::Vec k, l, m;
  _lineStream >> k[0] >> l[0] >> m[0]
          >> k[1] >> l[1] >> m[1]
          >> k[2] >> l[2] >> m[2];
  sys.box().K() = k;
  sys.box().L() = l;
  sys.box().M() = m;

  sys.box().update_triclinic();

  if (_lineStream.fail())
    throw InG96::Exception("Bad line in TRICLINICBOX block:\n" + *it +
          "\nTrying to read nine doubles");

}

void gio::InG96_i::readGenbox(System &sys) {
  std::vector<std::string> buffer;

  getblock(buffer);
  if (buffer[buffer.size() - 1].find("END") != 0)
    throw InG96::Exception("Coordinate file " + name() +
          " is corrupted. No END in GENBOX"
          " block. Got\n"
          + buffer[buffer.size() - 1]);

  std::string s;
  gio::concatenate(buffer.begin(), buffer.end() - 1, s);

  _lineStream.clear();
  _lineStream.str(s);
  int ntb;
  _lineStream >> ntb;
  if (_lineStream.fail())
    throw InG96::Exception("Bad line in GENBOX block:\n" + s +
          "\nTrying to read NTB integer.");
  double a, b, c, alpha, beta, gamma, phi, theta, psi, X, Y, Z;
  _lineStream >> a >> b >> c
          >> alpha >> beta >> gamma
          >> phi >> theta >> psi;

  if (_lineStream.fail())
    throw InG96::Exception("Bad line in GENBOX block:\n" + s +
          "\nTrying to read twelve doubles");
  _lineStream >> X >> Y >> Z;
  if (_lineStream.fail())
    X = Y = Z = 0.0;

  // Gromos++ is not implemented for X = Y = Z != 0
  if (X != 0.0 || Y != 0.0 || Z != 0.0) {
    throw InG96::Exception("GROMOS++ is not implemented for X = Y = Z != 0");
  }

  sys.box() = Box(gcore::Box::boxshape_enum(ntb),
          a, b, c, alpha, beta, gamma, phi, theta, psi, X, Y, Z);
}

void gio::InG96_i::readRemd(gcore::System &sys) {
  std::vector<std::string> buffer;
  getblock(buffer);
  if (buffer[buffer.size() - 1].find("END") != 0)
    throw InG96::Exception("Coordinate file " + name() +
          " is corrupted. No END in REMD"
          " block. Got\n"
          + buffer[buffer.size() - 1]);

  std::string s;

  gio::concatenate(buffer.begin(), buffer.end() - 1, s);

  _lineStream.clear();
  _lineStream.str(s);
  _lineStream >> sys.remd().id()
          >> sys.remd().run()
          >> sys.remd().temperature()
          >> sys.remd().lambda()
          >> sys.remd().Ti()
          >> sys.remd().li()
          >> sys.remd().Tj()
          >> sys.remd().lj()
          >> sys.remd().reeval();
  if (_lineStream.fail())
    throw InG96::Exception("Bad line in REMD block:\n" + s +
          "\nTrying to read all data");
}

InG96 &InG96::operator>>(System &sys) {

  if (!d_this) {
    throw InG96::Exception("read in frame, but no open file");
  }

  if (!d_this->stream())
    throw Exception("File " + name() + " is corrupted.");

  const std::string first = d_this->d_current;
  // std::cerr << first << std::endl;
  std::vector<std::string> buffer;
  bool readpos = false;
  bool readlatticeshifts = false;
  bool readvel = false;
  bool readcosDisplacement = false;
  bool readbox = false;
  bool readremd = false;
  
  d_this->initializeBfactors(sys);

  // skip frames
  // std::cerr << "operator<< : skip=" << d_this->d_skip << std::endl;

  for (; d_this->d_skip > 0; --d_this->d_skip) {

    do {
      // std::cerr << "skipping block " << d_this->d_current << std::endl;

      d_this->skipblock();
      d_this->getline(d_this->d_current);

    } while (d_this->d_current != first &&
            (!d_this->stream().eof()));

    if (d_this->stream().eof()) {
      // std::cerr << "skip eof: " << d_this->d_skip << std::endl;

      --d_this->d_skip;
      d_stride_eof = true;
      return *this;
    }

  }

  // only stride if not skip because of eof during last stride
  if (d_stride_eof == false) {
    int i = 1;
    for (; i < d_this->d_stride; ++i) {

      do {

        d_this->skipblock();
        d_this->getline(d_this->d_current);

      } while (d_this->d_current != first &&
              (!d_this->stream().eof()));

      if (d_this->stream().eof()) {
        // save remaining strides in skip for next file

        std::cerr << "stride eof: " << d_this->d_stride
                << "\ti: " << i << std::endl;

        d_this->d_skip = d_this->d_stride - i - 1;
        d_stride_eof = true;
        return *this;
      }

    }
  }

  // skipping and striding worked...
  d_stride_eof = false;

  do {
    switch (BLOCKTYPE[d_this->d_current]) {
      case titleblock:
        break;
      case timestep:
        d_this->readTimestep();
        break;
      case positionred:
        d_this->readPosition(sys);
        readpos = true;
        break;
      case position:
        d_this->readPosition(sys);
        readpos = true;
        break;
      case latticeshifts:
        d_this->readLatticeshifts(sys);
        readlatticeshifts = true;
        break;
      case velocityred:
        d_this->readVelocity(sys);
        readvel = true;
        break;
      case velocity:
        d_this->readVelocity(sys);
        readvel = true;
        break;
      case cosdisplacements:
        d_this->readCosDisplacements(sys);
        readcosDisplacement = true;
        break;
      case box:
        d_this->readBox(sys);
        readbox = true;
        break;
      case triclinicbox:
        d_this->readTriclinicbox(sys);
        readbox = true;
        break;
      case gmxbox:
        d_this->readGmxbox(sys);
        readbox = true;
        break;
      case genbox:
        d_this->readGenbox(sys);
        readbox = true;
        break;
      case remd:
        d_this->readRemd(sys);
        readremd = true;
        break;
      default:
        throw
        Exception("Block " + d_this->d_current +
                " is unknown in a coordinate file");
        break;
    }
    d_this->getline(d_this->d_current);
  } while (d_this->d_current != first && !d_this->stream().eof());

  sys.hasPos = readpos;
  sys.hasLatticeshifts = readlatticeshifts;
  sys.hasVel = readvel;
  sys.hasCosDisplacements = readcosDisplacement;
  sys.hasBox = readbox;
  sys.hasRemd = readremd;
  return *this;
}

InG96 &InG96::operator>>(utils::Time &time) {
  if (time.read()) {
    // we read the time from the traj: so just set it
    if (!d_this->d_time_read) {
      throw Exception("trajectory does not contain a TIMESTEP block. "
              "Use @time. ");
    }
    time.time() = d_this->d_time;
    time.steps() = d_this->d_step;
  } else {
    // we have to calculate the time
    time.time() += time.dt();
    time.steps()++; 
  }

  return *this;
}

std::string InG96::name()const {
  return d_this->name();
}

int InG96::skip()const {
  if (d_this)
    return d_this->d_skip;
  else return d_skip;
}

int InG96::stride()const {
  if (d_this)
    return d_this->d_stride;
  else return d_stride;
}

void InG96::skip(int skip) {
  if (d_this)
    d_this->d_skip = skip;
  else
    d_skip = skip;
}

void InG96::stride(int stride) {
  if (d_this)
    d_this->d_stride = stride;
  else
    d_stride = stride;
}

bool InG96::stride_eof()const {
  return d_stride_eof;
}
