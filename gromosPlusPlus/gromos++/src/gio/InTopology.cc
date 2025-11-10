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

// gio_InTopology.cc
#include "InTopology.h"

#include <cassert>
#include <iostream>
#include <map>
#include <map>
#include <cmath>
#include <iomanip>
#include <math.h>
#include <sstream>
#include <string>
#include <vector>

#include "Ginstream.h"
#include "../gromos/Exception.h"
#include "../gcore/VirtualAtomType.h"
#include "../gcore/BondType.h"
#include "../gcore/Bond.h"
#include "../gcore/AngleType.h"
#include "../gcore/Angle.h"
#include "../gcore/Constraint.h"
#include "../gcore/DihedralType.h"
#include "../gcore/Dihedral.h"
#include "../gcore/CrossDihedral.h"
#include "../gcore/ImproperType.h"
#include "../gcore/Improper.h"
#include "../gcore/LJException.h"
#include "../gcore/LJExceptionType.h"
#include "../gcore/LJType.h"
#include "../gcore/CGType.h"
#include "../gcore/AtomPair.h"
#include "../gcore/Exclusion.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Molecule.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/System.h"
#include "../gcore/GromosForceField.h"
#include "../gcore/LinearTopology.h"
#include "../args/Arguments.h"
#include "../gmath/Physics.h"

using namespace std;
using namespace gcore;
using gio::InTopology_i;
using gio::InTopology;

// Implementation class

class gio::InTopology_i : public gio::Ginstream {
    friend class InTopology;
    gcore::GromosForceField d_gff;
    gcore::System d_sys;
    std::string d_version;
    std::map<std::string, std::vector<std::string> > d_blocks;
    /**
     * the init function reads in the whole file into the map of blocks and
     * reads in the topology version
     */
    void init();
    /**
     * parseForceField takes all blocks that end up in the forcefield
     * and stores the information in... d_gff
     */
    void parseForceField();
    /**
     * parseSystem takes all blocks that end up in the system
     * and stores the information in d_sys;
     */
    void parseSystem();
    /**
     * _initBlock is a function that initialized the reading of
     * a block. It checks whether the block is read in, and returns
     * the number of elements that are to be read in from it.
     */
    int _initBlock(std::vector<std::string> &buffer,
            std::vector<std::string>::const_iterator &it,
            const std::string blockname,
            bool required);

    InTopology_i(std::string &s) : d_gff(), d_sys(), d_version(), d_blocks() {
        this->open(s);
    }

    ~InTopology_i(){
        this->close();
    }
};

gio::InTopology::InTopology(std::string name) {
    d_this = new InTopology_i(name);
    d_this->init();
    d_this->parseForceField();
    d_this->parseSystem();
}

gio::InTopology::~InTopology() {
    delete d_this;
}

const string &InTopology::version()const {
    return d_this->d_version;
}

const string InTopology::title()const {
    return d_this->title();
}

const gcore::System &InTopology::system()const {
    return d_this->d_sys;
}

const GromosForceField &InTopology::forceField()const {
    return d_this->d_gff;
}

int gio::InTopology_i::_initBlock(std::vector<std::string> &buffer,
        std::vector<std::string>::const_iterator &it,
        const string blockname,
        bool required=true) {
    int n;

    buffer.clear();
    buffer = d_blocks[blockname];
    if (buffer.size() < 3) {
        if (required){
          throw InTopology::Exception("Topology file " + name() +
            " is corrupted. No (or empty) " + blockname +
            " block!");
        } else {
          return 0;
        }
    }

    if (buffer[buffer.size() - 1].find("END") != 0)
        throw InTopology::Exception("Topology file " + name() +
            " is corrupted. No END in " + blockname +
            " block. Got\n"
            + buffer[buffer.size() - 1]);

    it = buffer.begin() + 1;
    _lineStream.clear();
    _lineStream.str(*it);
    _lineStream >> n;
    if (_lineStream.fail() || !_lineStream.eof()) {
        std::ostringstream msg;
        msg << "In block " << blockname << ": could not read number of lines. This "
                "number has to be the only value in the first line of the block.";
        throw gromos::Exception("InTopology", msg.str());
    }

    ++it;
    return n;
}

void gio::InTopology_i::init() {

    if (!stream())
        throw InTopology::Exception("Could not open topology file." + name());

    // First read the whole file into the map
    std::vector<std::string> buffer;
    std::vector<std::string>::const_iterator it;

    while (!stream().eof()) {
        getblock(buffer);
        if (buffer.size()) {
            d_blocks[buffer[0]] = buffer;
            buffer.clear();
        }
    }
    { // Version Block
        buffer.clear();
        buffer = d_blocks["TOPVERSION"];

        if (buffer.size() < 3)
            throw InTopology::Exception("Topology file " + name() +
                " is corrupted. No (or empty) TOPVERSION block!");
        if (buffer[buffer.size() - 1].find("END") != 0)
            throw InTopology::Exception("Topology file " + name() +
                " is corrupted. No END in TOPVERSION"
                " block. Got\n"
                + buffer[buffer.size() - 1]);
        it = buffer.begin() + 1;
        _lineStream.clear();
        _lineStream.str(*it);

        _lineStream >> d_version;
        if (_lineStream.fail())
            throw InTopology::Exception("Bad line in TOPVERSION block:\n" + *it);
    } // TOPVERSION
}

void gio::InTopology_i::parseForceField() {
    // generic variables
    double d[4];
    int i[5], num, n;
    string s;
    std::vector<std::string> buffer;
    std::vector<std::string>::const_iterator it;

    if (args::Arguments::inG96 == true) { // Topphyscon block:
        buffer.clear();
        buffer = d_blocks["TOPPHYSCON"];
        if (buffer.size() < 3)
            throw InTopology::Exception("Topology file " + name() +
                " is corrupted. No (or empty) TOPPHYSCON block!");
        if (buffer[buffer.size() - 1].find("END") != 0)
            throw InTopology::Exception("Topology file " + name() +
                " is corrupted. No END in TOPHYSCON"
                " block. Got\n"
                + buffer[buffer.size() - 1]);

        // this block comes as two lines with one number or as one line with two
        std::string topphyscon;
        gio::concatenate(buffer.begin() + 1, buffer.end() - 1, topphyscon);
        _lineStream.clear();
        _lineStream.str(topphyscon);
        _lineStream >> d[0] >> d[1];

        if (_lineStream.fail())
            throw InTopology::Exception("Bad line in TOPPHYSCON block:\n" + topphyscon);
        gmath::physConst.set_four_pi_eps_i(d[0]);
        gmath::physConst.set_hbar(d[1]);
        gmath::physConst.calc();
        d_gff.setFpepsi(gmath::physConst.get_four_pi_eps_i());
        d_gff.setHbar(gmath::physConst.get_hbar());
        d_gff.setSpdl(gmath::physConst.get_speed_of_light());
        d_gff.setBoltz(gmath::physConst.get_boltzmann());
    }// TOPPHYSCON
    else { // PHYSICALCONSTANTS block
        buffer.clear();
        buffer = d_blocks["PHYSICALCONSTANTS"];
        if (buffer.size() < 3)
            throw InTopology::Exception("Topology file " + name() +
                " is corrupted. No (or empty) PHYSICALCONSTANTS block!");
        if (buffer[buffer.size() - 1].find("END") != 0)
            throw InTopology::Exception("Topology file " + name() +
                " is corrupted. No END in PHYSICALCONSTANTS"
                " block. Got\n"
                + buffer[buffer.size() - 1]);

        // this block comes as four lines with one number or as one line with four or ...
        std::string physicalconstants;
        gio::concatenate(buffer.begin() + 1, buffer.end() - 1, physicalconstants);
        _lineStream.clear();
        _lineStream.str(physicalconstants);
        _lineStream >> d[0] >> d[1] >> d[2] >> d[3];

        if (_lineStream.fail())
            throw InTopology::Exception("Bad line in PHYSICALCONSTANTS block:\n" + physicalconstants);
        gmath::physConst.set_four_pi_eps_i(d[0]);
        gmath::physConst.set_hbar(d[1]);
        gmath::physConst.set_speed_of_light(d[2]);
        gmath::physConst.set_boltzmann(d[3]);
        gmath::physConst.calc();
        d_gff.setFpepsi(gmath::physConst.get_four_pi_eps_i());
        d_gff.setHbar(gmath::physConst.get_hbar());
        d_gff.setSpdl(gmath::physConst.get_speed_of_light());
        d_gff.setBoltz(gmath::physConst.get_boltzmann());
  } // PHYSICALCONSTANTS

    { // AtomTypename
        num = _initBlock(buffer, it, "ATOMTYPENAME");
        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
            _lineStream.clear();
            _lineStream.str(*it);
            _lineStream >> s;
            if (_lineStream.fail())
                throw InTopology::Exception("Bad line in ATOMTYPENAME block:\n" + *it);
            d_gff.addAtomTypeName(s);
        }
        if (n != num) {
            ostringstream os;
            os << "Incorrect number of AtomTypes in ATOMTYPENAME block\n"
                    << "Expected " << num << ", but found " << n;
            throw InTopology::Exception(os.str());
        }

    } // ATOMTYPENAME
    if (d_blocks["VIRTUALATOMTYPE"].size() > 2) {
 
    // temporary vectors for force constants and bond lengths

      buffer.clear();
      num = _initBlock(buffer, it, "VIRTUALATOMTYPE");
      for (n = 0; it != buffer.end() - 1; ++it, ++n) {
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> i[0] >> d[0] >> d[1];
        if (_lineStream.fail())
           throw InTopology::Exception("Bad line in VIRTUALATOMTYPE block:\n" + *it);
        d_gff.addVirtualAtomType(VirtualAtomType(i[0], d[0], d[1]));
      }
      if (n != num) {
        ostringstream os;
        os << "Incorrect number of Virtual Atom Types in VIRTUALATOMTYPE block\n"
           << "Expected " << num << ", but found " << n;
        throw InTopology::Exception(os.str());
      }
    } // VIRTUALATOMTYPE
 
    { // bond types
    // BONDTYPE, HARMBONDTYPE and BONDSTRETCHTYPE blocks

    // temporary vectors for force constants and bond lengths
    std::vector<double> harm_k, harm_d0, quart_k, quart_d0;
    int foundBondtypeBlock = 0;
    int foundBondstretchtypeBlock = 0;
    int foundHarmbondtypeBlock = 0;
    if (args::Arguments::inG96 == false) {

        buffer.clear();
        buffer = d_blocks["BONDTYPE"];
        if (buffer.size() >= 3) foundBondtypeBlock = 1;

        buffer.clear();
        buffer = d_blocks["BONDSTRETCHTYPE"];
        if (buffer.size() >= 3) foundBondstretchtypeBlock = 1;

        buffer.clear();
        buffer = d_blocks["HARMBONDTYPE"];
        if (buffer.size() >= 3) foundHarmbondtypeBlock = 1;

        // 
        if (foundBondtypeBlock == 0 && foundBondstretchtypeBlock == 0 && foundHarmbondtypeBlock == 0) {
            throw InTopology::Exception("Corrupt topology file:\n"
                    "all of BONDTYPE, HARMBONDTYPE and BONDSTRETCHTYPE block are missing or empty!");
        }
        
        // do not allow BONDTYPE/HARMBONDTYPE when there is a BONDSTRETCHTYPE block
        if ((foundBondtypeBlock == 1 || foundHarmbondtypeBlock == 1) && foundBondstretchtypeBlock == 1) {
            throw InTopology::Exception("Topology file:\n"
                    "use either (BONDTYPE and/or HARMBONDTYPE) OR BONDSTRETCHTYPE block, but do not mix them!");
        }
    }

    if (foundBondtypeBlock) {
        std::cerr << "# reading bond types from BONDTYPE block\n";
        num = _initBlock(buffer, it, "BONDTYPE");
        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
            _lineStream.clear();
            _lineStream.str(*it);
            _lineStream >> d[0] >> d[1];
            if (_lineStream.fail())
                throw InTopology::Exception("Bad line in BONDTYPE block:\n" + *it);

            quart_k.push_back(d[0]);
            quart_d0.push_back(d[1]);
            //d_gff.addBondType(BondType(n, d[0], d[1]));
        }
        if (n != num) {
            ostringstream os;
            os << "Incorrect number of BondTypes in BONDTYPE block\n"
                    << "Expected " << num << ", but found " << n;
            throw InTopology::Exception(os.str());
        }
    } // BONDTYPE
    
    if (foundHarmbondtypeBlock) {
        std::cerr << "# reading bond types from HARMBONDTYPE block\n";
        num = _initBlock(buffer, it, "HARMBONDTYPE");
        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
            _lineStream.clear();
            _lineStream.str(*it);
            _lineStream >> d[0] >> d[1];
            if (_lineStream.fail())
                throw InTopology::Exception("Bad line in HARMBONDTYPE block:\n" + *it);

            harm_k.push_back(d[0]);
            harm_d0.push_back(d[1]);

            //d_gff.addBondType(BondType(n, d[0], d[1],false));
        }
        if (n != num) {
            ostringstream os;
            os << "Incorrect number of BondTypes in HARMBONDTYPE block\n"
                    << "Expected " << num << ", but found " << n;
            throw InTopology::Exception(os.str());
        }
    } // HARMBONDTYPE
    
    if (foundBondstretchtypeBlock) { // BONDSTRETCHTYPE block
        num = _initBlock(buffer, it, "BONDSTRETCHTYPE");
        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
            _lineStream.clear();
            _lineStream.str(*it);
            _lineStream >> d[0] >> d[1] >> d[2];
            if (_lineStream.fail())
                throw InTopology::Exception("Bad line in BONDSTRETCHTYPE block:\n" + *it);
            //d_gff.addBondType(BondType(n, d[0], d[1], d[2]));
            harm_k.push_back(d[1]);
            harm_d0.push_back(d[2]);
            quart_k.push_back(d[0]);
            quart_d0.push_back(d[2]);
        }
        if (n != num) {
            ostringstream os;
            os << "Incorrect number of BondTypes in BONDSTRETCHTYPE block\n"
                    << "Expected " << num << ", but found " << n;
            throw InTopology::Exception(os.str());
        }
    } // BONDSTRETCHTYPE
    
    if (quart_k.size() != 0 && harm_k.size() != 0 && quart_k.size() != harm_k.size()) {
          throw InTopology::Exception("Different number of bond types in BONDTYPE and HARMBONDTYPE block.");
    }
    if (quart_k.size() && harm_k.size()){
        for (unsigned int i=0; i < quart_k.size(); i++) {  
          if ( quart_d0[i] !=  harm_d0[i] ) {
            ostringstream msg;
            msg << "Bond type "<< i <<" in BONDTYPE (b0="<<quart_d0[i]<<") and HARMBONDTYPE (b0="<<harm_d0[i]<<") block do not agree.";
            throw InTopology::Exception(msg.str());
          }
          // no check if the harmonic and quartic force constants correspond to each other
          d_gff.addBondType(BondType(i, quart_k[i], harm_k[i], quart_d0[i]));
        }
    } else if (quart_k.size()) {  
          for (unsigned int i=0; i < quart_k.size(); i++) { 
            d_gff.addBondType(BondType(i, quart_k[i], quart_d0[i], true));
          }
    } else if (harm_k.size()) {
          for (unsigned int i=0; i < harm_k.size(); i++) { 
            d_gff.addBondType(BondType(i, harm_k[i], harm_d0[i], false));
          }
    }
    } // bond types

    
    { // BONDANGLETYPE, HARMBONDANGLETYPE and BONDANGLEBENDTYPE blocks
      
    // temporary vectors for force constants and angles
    std::vector<double> harm_k, harm_th, cosharm_k, cosharm_th;
    
    int foundBondangletypeBlock = 0;
    int foundHarmBondangletypeBlock = 0;
    int foundBondanglebendtypeBlock = 0;
    if (args::Arguments::inG96 == false) {

        buffer.clear();
        buffer = d_blocks["BONDANGLETYPE"];
        if (buffer.size() >= 3) foundBondangletypeBlock = 1;

        buffer.clear();
        buffer = d_blocks["BONDANGLEBENDTYPE"];
        if (buffer.size() >= 3) foundBondanglebendtypeBlock = 1;

        buffer.clear();
        buffer = d_blocks["HARMBONDANGLETYPE"];
        if (buffer.size() >= 3) foundHarmBondangletypeBlock = 1;

        // The next if statement is necessary so that we can
        // skip reading in of BONDANGLETYPE block if foundBondanglebendtypeBlock = 1
        if (foundBondangletypeBlock == 0 && foundBondanglebendtypeBlock == 0 && foundHarmBondangletypeBlock == 0) {
            throw InTopology::Exception("Corrupt topology file:\n"
                    "all of BONDANGLETYPE, HARMBONDANGLETYPE and BONDANGLEBENDTYPE block are missing or empty!");
        }
        if ((foundBondangletypeBlock == 1 || foundHarmBondangletypeBlock == 1) && foundBondanglebendtypeBlock == 1) {
            throw InTopology::Exception("use either (BONDANGLETYPE and/or HARMBONDANGLETYPE) OR BONDANGLEBENDTYPE block, but do not mix them!");
        }
    }

    if (foundBondangletypeBlock) { // BONDANGLETYPE block
        num = _initBlock(buffer, it, "BONDANGLETYPE");
        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
            _lineStream.clear();
            _lineStream.str(*it);
            _lineStream >> d[0] >> d[1];
            cosharm_k.push_back(d[0]);
            cosharm_th.push_back(d[1]);
            if (_lineStream.fail())
                throw InTopology::Exception("Bad line in BONDANGLETYPE block:\n" + *it);
        }
        if (n != num) {
            ostringstream os;
            os << "Incorrect number of AngleTypes in BONDANGLETYPE block\n"
                    << "Expected " << num << ", but found " << n;
            throw InTopology::Exception(os.str());
        }
    }// BONDANGLETYPE
    
    if (foundHarmBondangletypeBlock) { // HARMBONDANGLETYPE block
        num = _initBlock(buffer, it, "HARMBONDANGLETYPE");
        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
            _lineStream.clear();
            _lineStream.str(*it);
            _lineStream >> d[0] >> d[1];
            harm_k.push_back(d[0]);
            harm_th.push_back(d[1]);
            if (_lineStream.fail())
                throw InTopology::Exception("Bad line in HARMBONDANGLETYPE block:\n" + *it);
        }
        if (n != num) {
            ostringstream os;
            os << "Incorrect number of AngleTypes in HARMBONDANGLETYPE block\n"
                    << "Expected " << num << ", but found " << n;
            throw InTopology::Exception(os.str());
        }
    }// HARMBONDANGLETYPE
    
    if (foundBondanglebendtypeBlock) { // BONDANGLEBENDTYPE block
        num = _initBlock(buffer, it, "BONDANGLEBENDTYPE");
        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
            _lineStream.clear();
            _lineStream.str(*it);
            _lineStream >> d[0] >> d[1] >> d[2];
            cosharm_k.push_back(d[0]);
            cosharm_th.push_back(d[2]);
            harm_k.push_back(d[1]);
            harm_th.push_back(d[2]);
            if (_lineStream.fail())
                throw InTopology::Exception("Bad line in BONDANGLEBENDTYPE block:\n" + *it);
        }
        if (n != num) {
            ostringstream os;
            os << "Incorrect number of AngleTypes in BONDANGLEBENDTYPE block\n"
                    << "Expected " << num << ", but found " << n;
            throw InTopology::Exception(os.str());
        }
    } // BONDANGLEBENDTYPE
    
      if (cosharm_k.size() != 0 && harm_k.size() != 0 && cosharm_k.size() != harm_k.size()) {
          throw InTopology::Exception("Different number of angle types in BONDANGLETYPE and HARMBONDANGLETYPE block.");
      }
      if (cosharm_k.size() && harm_k.size()){
        for (unsigned int i=0; i < cosharm_k.size(); i++) {  
          if ( cosharm_th[i] !=  harm_th[i] ) {
            ostringstream msg;
            msg << "Bond angle "<< i <<" in BONDANGLETYPE (theta="<<cosharm_th[i]<<") and HARMBONDANGLETYPE (theta="<<harm_th[i]<<") block do not agree.";
            throw InTopology::Exception(msg.str());
          }
          // no check if the harmonic and cosine harmonic force constants correspond to each other
          d_gff.addAngleType(AngleType(i, cosharm_k[i], harm_k[i], cosharm_th[i]));
        }
      } else if (cosharm_k.size()) {  
          for (unsigned int i=0; i < cosharm_k.size(); i++) {          
            try {
                d_gff.addAngleType(AngleType(i, cosharm_k[i], cosharm_th[i]));
            } catch (gromos::Exception & exp) {
                if (!args::Arguments::outG96) {
                    std::cerr << exp.what() << std::endl
                            << "Setting harmonic force constant to -1.0." << std::endl;
                }
                d_gff.addAngleType(AngleType(i, cosharm_k[i], -1, cosharm_th[i]));
            }
          }
      } else if (harm_k.size()) {
            std::cerr << "# Warning: Angle types: only harmonic force constants (HARMBONDANGLETYPE block) found, setting cosine harmonic force constant to -1.0." << std::endl;
          for (unsigned int i=0; i < harm_k.size(); i++) { 
            d_gff.addAngleType(AngleType(i, -1, harm_k[i], harm_th[i]));
          }
      }
    } // bond angle types

    { // IMPDIHEDRALTYPE block
        num = _initBlock(buffer, it, "IMPDIHEDRALTYPE");
        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
            _lineStream.clear();
            _lineStream.str(*it);
            _lineStream >> d[0] >> d[1];
            if (_lineStream.fail())
                throw InTopology::Exception("Bad line in IMPDIHEDRALTYPE block:\n" + *it);
            d_gff.addImproperType(ImproperType(n, d[0], d[1]));
        }
        if (n != num) {
            ostringstream os;
            os << "Incorrect number of ImproperTypes in IMPDIHEDRALTYPE block\n"
                    << "Expected " << num << ", but found " << n;
            throw InTopology::Exception(os.str());
        }
    } // IMPDIHEDRALTYPE

    // DIHEDRALTYPE and TORSDIHEDRALTYPE blocks
    // GROMOS08: check which of the two blocks is given
    int foundDihedraltypeBlock = 0;
    int foundTorsdihedraltypeBlock = 0;
    if (args::Arguments::inG96 == false) {

        buffer.clear();
        buffer = d_blocks["DIHEDRALTYPE"];
        if (buffer.size() >= 3) foundDihedraltypeBlock = 1;

        buffer.clear();
        buffer = d_blocks["TORSDIHEDRALTYPE"];
        if (buffer.size() >= 3) foundTorsdihedraltypeBlock = 1;

        // The next if statement is necessary so that we can
        // skip reading in of DIHEDRALTYPE block if foundTorsdihedraltypeBlock = 1
        if (foundDihedraltypeBlock == 0 && foundTorsdihedraltypeBlock == 0) {
            throw InTopology::Exception("Corrupt topology file:\n"
                    "both DIHEDRALTYPE and TORSDIHEDRALTYPE block are missing or empty!");
        }
    }

    if (args::Arguments::inG96 == true || foundTorsdihedraltypeBlock == 0) { // DIHEDRALTYPE block
        num = _initBlock(buffer, it, "DIHEDRALTYPE");
        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
            _lineStream.clear();
            _lineStream.str(*it);
            _lineStream >> d[0] >> d[1] >> i[0];
            if (_lineStream.fail())
                throw InTopology::Exception("Bad line in DIHEDRALTYPE block:\n" + *it);
            d[2] = acos(d[1])*180.0 / M_PI;
            d_gff.addDihedralType(DihedralType(n, d[0], d[1], d[2], i[0]));
        }
        if (n != num) {
            ostringstream os;
            os << "Incorrect number of DihedralTypes in DIHEDRALTYPE block\n"
                    << "Expected " << num << ", but found " << n;
            throw InTopology::Exception(os.str());
        }
    }// DIHEDRALTYPE
  else { // TORSDIHEDRALTYPE block
    num = _initBlock(buffer, it, "TORSDIHEDRALTYPE");
    for (n = 0; it != buffer.end() - 1; ++it, ++n) {
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> d[0] >> d[2] >> i[0];
      if (_lineStream.fail())
        throw InTopology::Exception("Bad line in TORSDIHEDRALTYPE block:\n" + *it);
      d[1] = cos(d[2] * M_PI / 180.0);
      d_gff.addDihedralType(DihedralType(n, d[0], d[1], d[2], i[0]));
    }
    if (n != num) {
      ostringstream os;
      os << "Incorrect number of DihedralTypes in DIHEDRALTYPE block\n"
              << "Expected " << num << ", but found " << n;
      throw InTopology::Exception(os.str());
    }
  } // TORSDIHEDRALTYPE

  { // LJPARAMETERS
    num = _initBlock(buffer, it, "LJPARAMETERS");
    for (n = 0; it < buffer.end() - 1; ++it, ++n) {
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i[0] >> i[1] >> d[0] >> d[1] >> d[2] >> d[3];
      if (_lineStream.fail())
        throw InTopology::Exception("Bad line in LJPARAMETERS block:\n" + *it);
      d_gff.setLJType(AtomPair(--i[0], --i[1]), LJType(d[0], d[1], d[2], d[3]));
    }
    if (n != num) {
      ostringstream os;
      os << "Incorrect number of LJ parameters in LJPARAMETERS block\n"
              << "Expected " << num << ", but found " << n;
      throw InTopology::Exception(os.str());
    }
  } // LJPARAMETERS

  if (args::Arguments::inG96 != true) { // LJEXCEPTIONS
    num = _initBlock(buffer, it, "LJEXCEPTIONS");
    for (n = 0; it < buffer.end() - 1; ++it, ++n) {
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i[0] >> i[1] >> d[0] >> d[1];
      if (_lineStream.fail())
        throw InTopology::Exception("Bad line in LJEXCEPTIONS block:\n" + *it);
      LJExceptionType lj(n,d[0],d[1]);
      // save the thing twice: once for writing the topology (so the function also
      // works when building the topology form a bb and parameter file) ...
      d_gff.addLJExceptionType(lj);
      // ... and once for easyer accsess of the LJ exceptions, e.g. for energy calculations
      d_gff.setLJException(AtomPair(--i[0], --i[1]), lj);
    }
    if (n != num) {
      ostringstream os;
      os << "Incorrect number of LJ exceptions in LJEXCEPTIONS block\n"
              << "Expected " << num << ", but found " << n;
      throw InTopology::Exception(os.str());
    }
  }

  { // CGPARAMETERS
    // first check if the block is present at all
    if (d_blocks["CGPARAMETERS"].size() > 2) {
      num = _initBlock(buffer, it, "CGPARAMETERS");
      for (n = 0; it < buffer.end() - 1; ++it, ++n) {
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> i[0] >> i[1] >> d[0] >> d[1];
        if (_lineStream.fail())
          throw InTopology::Exception("Bad line in CGPARAMETERS block:\n" + *it);
        d_gff.setCGType(AtomPair(--i[0], --i[1]), CGType(d[0], d[1]));
      }
      if (n != num) {
        ostringstream os;
        os << "Incorrect number of coarse grain LJ parameters in CGPARAMETERS block\n"
                << "Expected " << num << ", but found " << n;
        throw InTopology::Exception(os.str());
      }
    }
  } // CGPARAMETERS
}

void gio::InTopology_i::parseSystem() {
  // something to distinguish H atoms from others.
  double Hmass = 1.008;

  // generic variables
  double d[5];
  int i[9], num, n;
  string s;
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  // classes to store the data temporarily
  gcore::LinearTopology lt;
  gcore::SolventTopology st;

  { // RESNAME block
    num = _initBlock(buffer, it, "RESNAME");
    for (n = 0; it != buffer.end() - 1; ++it, ++n) {
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> s;
      if (_lineStream.fail())
        throw InTopology::Exception("Bad line in RESNAME block:\n" + *it);
      lt.setResName(n, s);
    }
    if (n != num) {
      ostringstream os;
      os << "Incorrect number of residues in RESNAME block\n"
              << "Expected " << num << ", but found " << n;
      throw InTopology::Exception(os.str());
    }
  } // RESNAME

  { // SOLUTEATOM block

    num = _initBlock(buffer, it, "SOLUTEATOM");
    // put the rest of the buffer into a single stream
    std::string soluteAtoms;
    std::vector<std::string>::const_iterator sAb = it, sAe = buffer.end() - 1;
    gio::concatenate(sAb, sAe, soluteAtoms);
    _lineStream.clear();
    _lineStream.str(soluteAtoms);

    for (n = 0; n < num; n++) {
      lt.addAtom(AtomTopology());
      _lineStream >> i[0];
      if (i[0] != n + 1)
        throw InTopology::Exception("Atom numbers are not sequential!");

      // residue number
      _lineStream >> i[0];
      lt.setResNum(n, --i[0]);
      // Atom Name
      _lineStream >> s;
      lt.atoms()[n].setName(s);
      // IAC
      _lineStream >> i[0];
      lt.atoms()[n].setIac(--i[0]);
      // mass
      _lineStream >> d[0];
      lt.atoms()[n].setMass(d[0]);
      if (d[0] == Hmass) lt.atoms()[n].setH(true);
      // charge
      _lineStream >> d[0];
      lt.atoms()[n].setCharge(d[0]);
      // charge group code
      _lineStream >> i[0];
      lt.atoms()[n].setChargeGroup(i[0]);
      // Exclusions 1-2 and 1-3
      Exclusion *e;
      e = new Exclusion();
      _lineStream >> i[0];
      for (int l = 0; l < i[0]; l++) {
        _lineStream >> i[1];
        e->insert(--i[1]);
      }
      lt.atoms()[n].setExclusion(*e);
      // Exclusions 1-4
      delete e;
      e = new Exclusion();
      _lineStream >> i[0];
      for (int l = 0; l < i[0]; l++) {
        _lineStream >> i[1];
        e->insert(--i[1]);
      }
      lt.atoms()[n].setExclusion14(*e);
      delete e;
      if (_lineStream.fail()) {
        ostringstream os;
        os << "Bad line in SOLUTEATOM block. Atom" << n + 1;
        throw InTopology::Exception(os.str());
      }
    }

  } // SOLUTEATOM
  { // SOLUTEPOLARISATION block
    if (d_blocks["SOLUTEPOLARISATION"].size() > 2) {
      num = _initBlock(buffer, it, "SOLUTEPOLARISATION");

      for (n = 0; it != buffer.end() - 1; ++it, ++n) {
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> i[0] >> d[0] >> d[1] >> d[2] >> d[3] >> d[4] >> i[1] >> i[2];

        if (_lineStream.fail()) {
          ostringstream os;
          os << "Bad line in SOLUTEPOLARISATION block. Line" << n + 1;
          throw InTopology::Exception(os.str());
        }
        i[0]--;
        i[1]--;
        i[2]--;

        if (i[0] < 0 || i[0] >= int(lt.atoms().size())) {
          ostringstream os;
          os << "SOLUTEPOLARISATION block: bad atom number " << i[0] + 1 << ".";
          throw InTopology::Exception(os.str());
        }
        if (i[1] < 0 || i[1] >= int(lt.atoms().size())) {
          ostringstream os;
          os << "SOLUTEPOLARISATION block: bad atom number " << i[1] + 1 << ".";
          throw InTopology::Exception(os.str());
        }
        if (i[2] < 0 || i[2] >= int(lt.atoms().size())) {
          ostringstream os;
          os << "SOLUTEPOLARISATION block: bad atom number " << i[2] + 1 << ".";
          throw InTopology::Exception(os.str());
        }

        lt.atoms()[i[0]].setPolarisable(true);
        lt.atoms()[i[0]].setPolarisability(d[0]);
        lt.atoms()[i[0]].setCosCharge(d[1]);
        lt.atoms()[i[0]].setDampingLevel(d[2]);
        lt.atoms()[i[0]].setDampingPower(d[3]);
        lt.atoms()[i[0]].setPoloffsiteGamma(d[4]);
        lt.atoms()[i[0]].setPoloffsiteI(i[1]);
        lt.atoms()[i[0]].setPoloffsiteJ(i[2]);
      } // for lines
      if (n != num) {
        ostringstream os;
        os << "Incorrect number of atoms in SOLUTEPOLARISATION block\n"
                << "Expected " << num << ", but found " << n;
        throw InTopology::Exception(os.str());
      }
    } // if block present
  } // SOLUTEPOLARISATION
  { // CGSOLUTE
    // first check if the block is present at all
    if (d_blocks["CGSOLUTE"].size() > 2) {
      num = _initBlock(buffer, it, "CGSOLUTE");
      for (n = 0; it < buffer.end() - 1; ++it, ++n) {
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> i[0] >> i[1] >> i[2];
        if (_lineStream.fail())
          throw InTopology::Exception("Bad line in CGSOLUTE block:\n" + *it);
        for (unsigned int j = (i[0]-1); int(j) < i[1]; ++j) {
          lt.atoms()[j].setCoarseGrained(true);
          lt.atoms()[j].setCGFactor(i[2]);
        }
      }
      if (n != num) {
        ostringstream os;
        os << "Incorrect number of ranges in CGSOLUTE block\n"
                << "Expected " << num << ", but found " << n;
        throw InTopology::Exception(os.str());
      }
    }
  } // CGSOLUTE
  { // BONDH
    num = _initBlock(buffer, it, "BONDH");
    for (n = 0; it < buffer.end() - 1; ++it, ++n) {
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i[0] >> i[1] >> i[2];
      if (_lineStream.fail())
        throw InTopology::Exception("Bad line in BONDH block:\n" + *it);
      Bond bond(--i[0], --i[1]);
      bond.setType(--i[2]);
      lt.addBond(bond);
    }
    if (n != num) {
      ostringstream os;
      os << "Incorrect number of bonds in BONDH block\n"
              << "Expected " << num << ", but found " << n;
      throw InTopology::Exception(os.str());
    }
  } // BONDH
  { // BOND
    num = _initBlock(buffer, it, "BOND");
    for (n = 0; it < buffer.end() - 1; ++it, ++n) {
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i[0] >> i[1] >> i[2];
      if (_lineStream.fail())
        throw InTopology::Exception("Bad line in BOND block:\n" + *it);
      Bond bond(--i[0], --i[1]);
      bond.setType(--i[2]);
      lt.addBond(bond);
    }
    if (n != num) {
      ostringstream os;
      os << "Incorrect number of bonds in BOND block\n"
              << "Expected " << num << ", but found " << n;
      throw InTopology::Exception(os.str());
    }
  } // BOND
  { // BONDDP
    // first check if the block is present at all
    if (d_blocks["BONDDP"].size() > 2) {
      num = _initBlock(buffer, it, "BONDDP");
      for (n = 0; it < buffer.end() - 1; ++it, ++n) {
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> i[0] >> i[1] >> i[2];
        if (_lineStream.fail())
          throw InTopology::Exception("Bad line in BONDDP block:\n" + *it);
        Bond bond(--i[0], --i[1]);
        bond.setType(--i[2]);
        lt.addDipoleBond(bond);
      }
      if (n != num) {
        ostringstream os;
        os << "Incorrect number of bonds in BONDDP block\n"
                << "Expected " << num << ", but found " << n;
        throw InTopology::Exception(os.str());
      }
    } // if block is present
  } // BONDDP

  { // BONDANGLEH
    num = _initBlock(buffer, it, "BONDANGLEH");
    for (n = 0; it < buffer.end() - 1; ++it, ++n) {
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i[0] >> i[1] >> i[2] >> i[3];
      if (_lineStream.fail())
        throw InTopology::Exception("Bad line in BONDANGLEH block:\n" + *it);
      Angle angle(--i[0], --i[1], --i[2]);
      angle.setType(--i[3]);
      lt.addAngle(angle);
    }
    if (n != num) {
      ostringstream os;
      os << "Incorrect number of angles in BONDANGLEH block\n"
              << "Expected " << num << ", but found " << n;
      throw InTopology::Exception(os.str());
    }
  } // BONDANGLEH
  { // BONDANGLE
    num = _initBlock(buffer, it, "BONDANGLE");
    for (n = 0; it < buffer.end() - 1; ++it, ++n) {
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i[0] >> i[1] >> i[2] >> i[3];
      if (_lineStream.fail())
        throw InTopology::Exception("Bad line in BONDANGLE block:\n" + *it);
      Angle angle(--i[0], --i[1], --i[2]);
      angle.setType(--i[3]);
      lt.addAngle(angle);
    }
    if (n != num) {
      ostringstream os;
      os << "Incorrect number of angles in BONDANGLE block\n"
              << "Expected " << num << ", but found " << n;
      throw InTopology::Exception(os.str());
    }
  } // BONDANGLE

  { // IMPDIHEDRALH
    num = _initBlock(buffer, it, "IMPDIHEDRALH");
    for (n = 0; it < buffer.end() - 1; ++it, ++n) {
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i[0] >> i[1] >> i[2] >> i[3] >> i[4];
      if (_lineStream.fail())
        throw InTopology::Exception("Bad line in IMPDIHEDRALH block:\n" + *it);
      Improper improper(--i[0], --i[1], --i[2], --i[3]);
      improper.setType(--i[4]);
      lt.addImproper(improper);
    }
    if (n != num) {
      ostringstream os;
      os << "Incorrect number of impropers in IMPDIHEDRALH block\n"
              << "Expected " << num << ", but found " << n;
      throw InTopology::Exception(os.str());
    }
  } // IMPDIHEDRALH
  { // IMPDIHEDRAL
    num = _initBlock(buffer, it, "IMPDIHEDRAL");
    for (n = 0; it < buffer.end() - 1; ++it, ++n) {
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i[0] >> i[1] >> i[2] >> i[3] >> i[4];
      if (_lineStream.fail())
        throw InTopology::Exception("Bad line in DIHEDRAL block:\n" + *it);
      Improper improper(--i[0], --i[1], --i[2], --i[3]);
      improper.setType(--i[4]);
      lt.addImproper(improper);
    }
    if (n != num) {
      ostringstream os;
      os << "Incorrect number of Impropers in IMPDIHEDRAL block\n"
              << "Expected " << num << ", but found " << n;
      throw InTopology::Exception(os.str());
    }
  } // IMPDIHEDRAL

  { // DIHEDRALH
    num = _initBlock(buffer, it, "DIHEDRALH");
    for (n = 0; it < buffer.end() - 1; ++it, ++n) {
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i[0] >> i[1] >> i[2] >> i[3] >> i[4];
      if (_lineStream.fail())
        throw InTopology::Exception("Bad line in DIHEDRALH block:\n" + *it);
      Dihedral dihedral(--i[0], --i[1], --i[2], --i[3]);
      dihedral.setType(--i[4]);
      lt.addDihedral(dihedral);
    }
    if (n != num) {
      ostringstream os;
      os << "Incorrect number of dihedrals in DIHEDRALH block\n"
              << "Expected " << num << ", but found " << n;
      throw InTopology::Exception(os.str());
    }
  } // DIHEDRALH
  { // DIHEDRAL
    num = _initBlock(buffer, it, "DIHEDRAL");
    for (n = 0; it < buffer.end() - 1; ++it, ++n) {
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i[0] >> i[1] >> i[2] >> i[3] >> i[4];
      if (_lineStream.fail())
        throw InTopology::Exception("Bad line in DIHEDRAL block:\n" + *it);
      Dihedral dihedral(--i[0], --i[1], --i[2], --i[3]);
      dihedral.setType(--i[4]);
      lt.addDihedral(dihedral);
    }
    if (n != num) {
      ostringstream os;
      os << "Incorrect number of dihedrals in DIHEDRAL block\n"
              << "Expected " << num << ", but found " << n;
      throw InTopology::Exception(os.str());
    }
  } // DIHEDRAL

  if (!args::Arguments::inG96) {
    if (d_blocks["CROSSDIHEDRALH"].size() > 2) { // CROSSDIHEDRALH
      num = _initBlock(buffer, it, "CROSSDIHEDRALH");
      for (n = 0; it < buffer.end() - 1; ++it, ++n) {
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> i[0] >> i[1] >> i[2] >> i[3] >> i[4] >> i[5] >> i[6] >> i[7] >> i[8];
        if (_lineStream.fail())
          throw InTopology::Exception("Bad line in CROSSDIHEDRALH block:\n" + *it);
        CrossDihedral crossdihedral(--i[0], --i[1], --i[2], --i[3], --i[4], --i[5], --i[6], --i[7]);
        crossdihedral.setType(--i[8]);
        lt.addCrossDihedral(crossdihedral);
      }
      if (n != num) {
        ostringstream os;
        os << "Incorrect number of cross dihedrals in CROSSDIHEDRALH block\n"
                << "Expected " << num << ", but found " << n;
        throw InTopology::Exception(os.str());
      }
    } // CROSSDIHEDRALH
    if (d_blocks["CROSSDIHEDRAL"].size() > 2) { // CROSSDIHEDRAL
      num = _initBlock(buffer, it, "CROSSDIHEDRAL");
      for (n = 0; it < buffer.end() - 1; ++it, ++n) {
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> i[0] >> i[1] >> i[2] >> i[3] >> i[4] >> i[5] >> i[6] >> i[7] >> i[8];
        if (_lineStream.fail())
          throw InTopology::Exception("Bad line in CROSSDIHEDRAL block:\n" + *it);
        CrossDihedral crossdihedral(--i[0], --i[1], --i[2], --i[3], --i[4], --i[5], --i[6], --i[7]);
        crossdihedral.setType(--i[8]);
        lt.addCrossDihedral(crossdihedral);
      }
      if (n != num) {
        ostringstream os;
        os << "Incorrect number of cross dihedrals in CROSSDIHEDRAL block\n"
                << "Expected " << num << ", but found " << n;
        throw InTopology::Exception(os.str());
      }
    } // CROSSDIHEDRAL
  }
  
  if (args::Arguments::inG96 != true){ // LJEXCEPTIONS
    num = _initBlock(buffer, it, "LJEXCEPTIONS");
    for (n = 0; it < buffer.end() - 1; ++it, ++n) {
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i[0] >> i[1] >> d[0] >> d[1];
      if (_lineStream.fail())
        throw InTopology::Exception("Bad line in LJEXCEPTIONS block:\n" + *it);
      LJException lj(--i[0], --i[1]);
      lj.setType(n);
      lt.addLJException(lj);
    }
    if (n != num) {
      ostringstream os;
      os << "Incorrect number of bonds in LJEXCEPTIONS block\n"
              << "Expected " << num << ", but found " << n;
      throw InTopology::Exception(os.str());
    }
  } // LJEXCEPTIONS

  { // SOLVENTATOM
    num = _initBlock(buffer, it, "SOLVENTATOM");
    if (num == 0)
      throw InTopology::Exception(
            "Cannot have 0 solvent atoms on the topology level!");
    std::string solventAtoms;
    std::vector<std::string>::const_iterator sAb = it, sAe = buffer.end() - 1;
    gio::concatenate(sAb, sAe, solventAtoms);
    _lineStream.clear();
    _lineStream.str(solventAtoms);
    for (n = 0; n < num; n++) {
      st.addAtom(AtomTopology());
      _lineStream >> i[0];
      if (i[0] != n + 1)
        throw InTopology::Exception(
              "Solvent atom numbers are not sequential!");
      // set name
      _lineStream >> s;
      st.atom(n).setName(s);
      // set IAC
      _lineStream >> i[0];
      st.atom(n).setIac(--i[0]);
      // set mass
      _lineStream >> d[0];
      st.atom(n).setMass(d[0]);
      // set charge
      _lineStream >> d[0];
      st.atom(n).setCharge(d[0]);

      if (_lineStream.fail()) {
        ostringstream os;
        os << "Bad line in SOLVENTATOM block. Atom " << n + 1;
        throw InTopology::Exception(os.str());
      }
    }
  } // SOLVENTATOM
  { // SOLVENTPOLARISATION block
    if (d_blocks["SOLVENTPOLARISATION"].size() > 2) {
      num = _initBlock(buffer, it, "SOLVENTPOLARISATION");
      for (n = 0; it < buffer.end() - 1; ++it, ++n) {
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> i[0] >> d[0] >> d[1] >> d[2] >> d[3] >> d[4] >> i[1] >> i[2];

        if (_lineStream.fail()) {
          ostringstream os;
          os << "Bad line in SOLVENTPOLARISATION block. Line" << n + 1;
          throw InTopology::Exception(os.str());
        }
        i[0]--;
        i[1]--;
        i[2]--;

        if (i[0] < 0 || i[0] >= st.numAtoms()) {
          ostringstream os;
          os << "SOLVENTPOLARISATION block: bad atom number " << i[0] + 1 << ".";
          throw InTopology::Exception(os.str());
        }
         if (i[1] < 0 || i[1] >= st.numAtoms()) {
          ostringstream os;
          os << "SOLVENTPOLARISATION block: bad atom number " << i[1] + 1 << ".";
          throw InTopology::Exception(os.str());
        }
         if (i[2] < 0 || i[2] >= st.numAtoms()) {
          ostringstream os;
          os << "SOLVENTPOLARISATION block: bad atom number " << i[2] + 1 << ".";
          throw InTopology::Exception(os.str());
        }

        st.atom(i[0]).setPolarisable(true);
        st.atom(i[0]).setPolarisability(d[0]);
        st.atom(i[0]).setCosCharge(d[1]);
        st.atom(i[0]).setDampingLevel(d[2]);
        st.atom(i[0]).setDampingPower(d[3]);
        st.atom(i[0]).setPoloffsiteGamma(d[4]);
        st.atom(i[0]).setPoloffsiteI(i[1]);
        st.atom(i[0]).setPoloffsiteJ(i[2]);

      } // for lines
      if (n != num) {
        ostringstream os;
        os << "Incorrect number of atoms in SOLVENTPOLARISATION block\n"
                << "Expected " << num << ", but found " << n;
        throw InTopology::Exception(os.str());
      }
    } // if block present
  } // SOLVENTPOLARISATION
  { // SOLVENTCONSTR
    num = _initBlock(buffer, it, "SOLVENTCONSTR");
    for (n = 0; it < buffer.end() - 1; ++it, ++n) {
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i[0] >> i[1] >> d[0];
      if (_lineStream.fail())
        throw InTopology::Exception("Bad line in SOLVENTCONSTR block:\n" + *it);
      Constraint constr(--i[0], --i[1]);
      constr.setDist(d[0]);
      st.addConstraint(constr);
    }
    if (n != num) {
      ostringstream os;
      os << "Incorrect number of constraints in SOLVENTCONSTR block\n"
              << "Expected " << num << ", but found " << n;
      throw InTopology::Exception(os.str());
    }
  } // SOLVENTCONSTR
  { // CONSTRAINT
    // this block is optional, check if it is there
    buffer.clear();
    buffer = d_blocks["CONSTRAINT"];
    num = _initBlock(buffer, it, "CONSTRAINT", false);
    for (n = 0; it < buffer.end() - 1; ++it, ++n) {
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i[0] >> i[1] >> i[2];
      if (_lineStream.fail())
        throw InTopology::Exception("Bad line in CONSTRAINT block:\n" + *it);
      Constraint constr(--i[0], --i[1]);
      constr.setType(--i[2]);
      double b0=d_gff.bondType(constr.bondtype()).b0();
      constr.setDist(b0);
      lt.addConstraint(constr); //TODO: check
    }
    if (n != num) {
      ostringstream os;
      os << "Incorrect number of constraints in CONSTRAINT block\n"
              << "Expected " << num << ", but found " << n;
      throw InTopology::Exception(os.str());
    }
  } // CONSTRAINT


  // Now parse the stuff into Topologies and the System.
  // this has to be done using this parse function of lt.parse
  // because we cannot use the System::operator= for a member function
  // (it involves a delete (this) statement
  lt.parse(d_sys);

  d_sys.addSolvent(Solvent(st));


  // virtual atoms rely on a system, so can only be parsed now
  if (d_blocks["VIRTUALATOMS"].size() > 2) {
    num = _initBlock(buffer, it, "VIRTUALATOMS");
    // put the rest of the buffer into a single stream
    std::string virtualAtoms;
    std::vector<std::string>::const_iterator sAb = it, sAe = buffer.end() - 1;
    gio::concatenate(sAb, sAe, virtualAtoms);
    _lineStream.clear();
    _lineStream.str(virtualAtoms);
    for (n = 0; n<num; ++n) {
      _lineStream >> i[0] >> i[1] >> d[0] >> i[2] >> i[3];
      std::vector<int> conf;
      for(int ii=0; ii< i[3]; ++ii){
        _lineStream >> i[4];
        conf.push_back(i[4]-1);
      }
      gcore::Exclusion e, e14;
      _lineStream >> i[3];
      for(int ii=0; ii <  i[3]; ++ii){
        _lineStream >> i[4];
        e.insert(i[4]-1);
      }
      _lineStream >> i[3];
      for(int ii=0; ii < i[3]; ++ii){
        _lineStream >> i[4];
        e14.insert(i[4]-1);
      }
      if (_lineStream.fail())
        throw InTopology::Exception("Bad line in VIRTUALATOMS block\n"+virtualAtoms);
      d_sys.addVirtualAtom(conf, i[2], 
                           d_gff.virtualAtomType(i[2]).dis1(), d_gff.virtualAtomType(i[2]).dis2(),
                           i[1]-1, d[0], e, e14);
    }
    if (n != num) {
      ostringstream os;
      os << "Incorrect number of virtual atoms in VIRTUALATOMS block\n"
              << "Expected " << num << ", but found " << n;
      throw InTopology::Exception(os.str());
    }
  }
  // In case of gromos08 topology, check if SOLUTEMOLECULE
  // block is consistent with connectivity. Crash if not.
  // Do this after parsing the system, it is only meant for checking...
  int totNumAt = 0;
  if (!(args::Arguments::inG96)) { // SOLUTEMOLECULES block
    num = _initBlock(buffer, it, "SOLUTEMOLECULES");
    _lineStream.clear();
    _lineStream.str(*it);

    if (num != d_sys.numMolecules()) {
      ostringstream os;
      os << "Incorrect number of solute molecules NSPM given in SOLUTEMOLECULES block\n"
              << "NSPM is set to " << i[0] << ", but from connectivity, I calculated "
              << d_sys.numMolecules() << " submolecules.\n\n"
              << "The block should look like this: \n";
      // SOLUTEMOLECULES block
      os << "SOLUTEMOLECULES\n"
              << "# NSPM: number of separate molecules in solute block\n"
              << "# NSP[1...NSPM]: atom sequence number of last atom\n"
              << "#                of the successive submolecules\n"
              << "#      NSPM  NSP[1...NSPM]\n";

      os << std::setw(10) << d_sys.numMolecules() << "\n";
      int nspmin = 0;
      for (int i = 0; i < d_sys.numMolecules(); ++i) {
        os << std::setw(6) << d_sys.mol(i).numAtoms() + nspmin;
        nspmin += d_sys.mol(i).numAtoms();
        if ((i + 1) % 10 == 0) os << "\n";
      }

      if (d_sys.numMolecules() % 10 != 0) os << "\n";
      os << "END\n";
      throw InTopology::Exception(os.str());
    }
    // put the rest of the buffer into a single stream
    std::string soluteMolecules;
    std::vector<std::string>::const_iterator sAb = it, sAe = buffer.end() - 1;
    gio::concatenate(sAb, sAe, soluteMolecules);
    _lineStream.clear();
    _lineStream.str(soluteMolecules);

    for (n = 0; n < num; n++) {
      totNumAt += d_sys.mol(n).numAtoms();
      _lineStream >> i[0];

      if (i[0] != totNumAt) {
        ostringstream os;
        os << "Incorrect number given in SOLUTEMOLECULES block for NSP[" << n + 1 << "]\n"
                << "got " << i[0] << ", but from connectivity, I calculated " << totNumAt << "\n\n"
                << "The block should look like this: \n";
        // SOLUTEMOLECULES block
        os << "SOLUTEMOLECULES\n"
                << "# NSPM: number of separate molecules in solute block\n"
                << "# NSP[1...NSPM]: atom sequence number of last atom\n"
                << "#                of the successive submolecules\n"
                << "#      NSPM  NSP[1...NSPM]\n";

        os << std::setw(10) << d_sys.numMolecules() << "\n";
        int nspmin = 0;
        for (int i = 0; i < d_sys.numMolecules(); ++i) {
          os << std::setw(6) << d_sys.mol(i).numAtoms() + nspmin;
          nspmin += d_sys.mol(i).numAtoms();
          if ((i + 1) % 10 == 0) os << "\n";
        }

        if (d_sys.numMolecules() % 10 != 0) os << "\n";
        os << "END\n";
        throw InTopology::Exception(os.str());
      }
    }
  } // SOLUTEMOLECULES


  // Add the temperature and pressure group blocks
  // (Set them to default values in case of gromos96
  if (args::Arguments::inG96 == true) {
    // how to get the total number of atoms in a tomology?!?
    // it is a bit odd but works at least, I hope there is something simpler
    // (even if I do not know where...)
    for (int m = 0; m < d_sys.numMolecules(); ++m) { // loop over molecules
      for (int a = 0; a < d_sys.mol(m).numAtoms(); ++a) { // loop over atoms in this molecule
        totNumAt++;
      }
    }
    d_sys.addTemperatureGroup(totNumAt);

    d_sys.addPressureGroup(totNumAt);
  } else {
    { // TEMPERATUREGROUPS block
      num = _initBlock(buffer, it, "TEMPERATUREGROUPS");
      _lineStream.clear();
      _lineStream.str(*it);
      // put the rest of the buffer into a single stream
      std::string temperatureGroups;
      std::vector<std::string>::const_iterator sAb = it, sAe = buffer.end() - 1;
      gio::concatenate(sAb, sAe, temperatureGroups);
      _lineStream.clear();
      _lineStream.str(temperatureGroups);

      for (n = 0; n < num; n++) {
        _lineStream >> i[0];
        if (_lineStream.fail())
          throw InTopology::Exception("Bad line in TEMPERATUREGROUPS block:\n" + *it);
        if (i[0] > totNumAt) {
          ostringstream os;
          os << "Error in TEMPERATUREGROUPS block:\n"
                  << "NST[" << n + 1 << "] > number of solute atoms!\n";
          throw InTopology::Exception(os.str());
        }
        if (n > 0) {
          if (i[0] == d_sys.temperatureGroup(n - 1)) {
            ostringstream os;
            os << "Last atom number in different TEMPERATUREGROUPS cannot be the same,\n"
                    << "but NST[" << n << "] = NST[" << n + 1 << "] !\n";
            throw InTopology::Exception(os.str());
          }
          if (i[0] < d_sys.temperatureGroup(n - 1)) {
            ostringstream os;
            os << "Values for NST[1...NSTM] in TEMPERATUREGROUPS should be sequential,\n"
                    << "but NST[" << n << "] > NST[" << n + 1 << "] !\n";
            throw InTopology::Exception(os.str());
          }
        }
        d_sys.addTemperatureGroup(i[0]);
      }
    } // TEMPERATUREGROUPS
    { // PRESSUREGROUPS block
      num = _initBlock(buffer, it, "PRESSUREGROUPS");
      _lineStream.clear();
      _lineStream.str(*it);
      // put the rest of the buffer into a single stream
      std::string pressureGroups;
      std::vector<std::string>::const_iterator sAb = it, sAe = buffer.end() - 1;
      gio::concatenate(sAb, sAe, pressureGroups);
      _lineStream.clear();
      _lineStream.str(pressureGroups);
      for (n = 0; n < num; ++n) {
        _lineStream >> i[0];
        if (_lineStream.fail())
          throw InTopology::Exception("Bad line in PRESSUREGROUPS block:\n" + *it);
        if (i[0] > totNumAt) {
          ostringstream os;
          os << "Error in PRESSUREGROUPS block:\n"
                  << "NSV[" << n + 1 << "] > number of solute atoms!\n";
          throw InTopology::Exception(os.str());
        }
        if (n > 0) {
          if (i[0] == d_sys.pressureGroup(n - 1)) {
            ostringstream os;
            os << "Last atom number in different PRESSUREGROUPS cannot be the same,\n"
                    << "but NSV[" << n << "] = NSV[" << n + 1 << "] !\n";
            throw InTopology::Exception(os.str());
          }
          if (i[0] < d_sys.pressureGroup(n - 1)) {
            ostringstream os;
            os << "Values for NSV[1...NSVM] in PRESSUREGROUPS should be sequential,\n"
                    << "but NSV[" << n << "] > NSV[" << n + 1 << "] !\n";
            throw InTopology::Exception(os.str());
          }
        }
        d_sys.addPressureGroup(i[0]);
      }
    } // PRESSUREGROUPS
  }

}

