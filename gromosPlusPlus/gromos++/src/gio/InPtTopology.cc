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
#include "InPtTopology.h"

#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <string>
#include <map>

#include "Ginstream.h"
#include "../gcore/Bond.h"
#include "../gcore/Angle.h"
#include "../gcore/Dihedral.h"
#include "../gcore/CrossDihedral.h"
#include "../gcore/Improper.h"
#include "../gcore/PtTopology.h"

using namespace std;
using namespace gcore;
using gio::InPtTopology_i;
using gio::InPtTopology;

// Implementation class
class gio::InPtTopology_i: public gio::Ginstream
{
  friend class InPtTopology;
  gcore::PtTopology d_pttopo;
  std::map<std::string, std::vector<std::string> > d_blocks;
  /**
   * the init function reads in the whole file into the map of blocks and 
   * reads in the topology version
   */
  void init();
  /**
   * parseForceField takes all blocks that end up in the forcefield
   * and stores the information in d_pttopo
   */
  void parsePtTopology();
  /**
   * _initBlock is a function that initialized the reading of
   * a block. It checks whether the block is read in, and returns 
   * the number of elements that are to be read in from it.
   */
  int _initBlock(std::vector<std::string> &buffer,
		 std::vector<std::string>::const_iterator &it,
		 const std::string blockname);

  InPtTopology_i(std::string &s):d_blocks()
  {
    this->open(s);
  }
  
};

gio::InPtTopology::InPtTopology(std::string name){
  d_this = new InPtTopology_i(name);
  d_this->init();
  d_this->parsePtTopology();
}

gio::InPtTopology::~InPtTopology()
{
  delete d_this;
}

const string InPtTopology::title()const
{
  return d_this->title();
}

const gcore::PtTopology &InPtTopology::ptTopo()const{
  return d_this->d_pttopo;
}

int gio::InPtTopology_i::_initBlock(std::vector<std::string> &buffer,
	       std::vector<std::string>::const_iterator &it,
	       const string blockname)
{
  int n;
  
  buffer.clear();
  buffer=d_blocks[blockname];
  if(buffer.size() < 3)
    throw InPtTopology::Exception("Topology file"+name()+
				" is corrupted. No (or empty) "+blockname+
				" block!");
  if(buffer[buffer.size()-1].find("END")!=0)
    throw InPtTopology::Exception("Topology file " + name() +
				     " is corrupted. No END in "+blockname+
				     " block. Got\n"
				     + buffer[buffer.size()-1]);

  it=buffer.begin() + 1;
  _lineStream.clear();
  _lineStream.str(*it);
  _lineStream >> n;
  ++it;
  return n;
}

void gio::InPtTopology_i::init(){
  
  if(!stream())
    throw InPtTopology::Exception("Could not open topology file."+name());
  
  // First read the whole file into the map
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  while(!stream().eof()){
    getblock(buffer);
    if(buffer.size()){
      d_blocks[buffer[0]]=buffer;
      buffer.clear();
    }
  }
}

void gio::InPtTopology_i::parsePtTopology() {
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  std::string nm;
  int k, l, n, iiac[2], numat, numpt;
  double mass[2], dq[2], alphaLJ, alphaCRF;

  if (d_blocks.count("PERTATOMPARAM")) {
    int numat = _initBlock(buffer, it, "PERTATOMPARAM");
    d_pttopo.setMultiPt(false);

    d_pttopo.setSize(numat, 2);
    d_pttopo.setPertName(0, "A");
    d_pttopo.setPertName(1, "B");
    for (n = 0; it < buffer.end() - 1; ++it, ++n) {
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> k >> l >> nm >> iiac[0] >> mass[0] >> dq[0]
              >> iiac[1] >> mass[1] >> dq[1] >> alphaLJ >> alphaCRF;

      if (_lineStream.fail()) {
        ostringstream os;
        os << "Bad line in PERTATOMPARAM block (line " << n + 1 << ").";
        throw InPtTopology::Exception(os.str());
      }

      d_pttopo.setAtomNum(n, k - 1);
      d_pttopo.setAtomName(n, nm);
      d_pttopo.setIac(n, 0, iiac[0] - 1);
      d_pttopo.setMass(n, 0, mass[0]);
      d_pttopo.setCharge(n, 0, dq[0]);
      d_pttopo.setIac(n, 1, iiac[1] - 1);
      d_pttopo.setMass(n, 1, mass[1]);
      d_pttopo.setCharge(n, 1, dq[1]);
      d_pttopo.setAlphaLJ(n, alphaLJ);
      d_pttopo.setAlphaCRF(n, alphaCRF);
    }
    if (n != numat)
      throw InPtTopology::Exception("Perturbation topology file " + name() +
            " is corrupted. Failed to read all atoms");

    if (d_blocks.count("PERTPOLPARAM")) {
      int numat = _initBlock(buffer, it, "PERTPOLPARAM");

      if (numat != d_pttopo.numAtoms()) {
        ostringstream os;
        os << "Error in PERTPOLPARAM block: The block contains " << numat
                << " atoms but the PERTATOMPARAM block " << d_pttopo.numAtoms()
                << " atoms.";
        throw InPtTopology::Exception(os.str());
      }
      for (n = 0; it < buffer.end() - 1; ++it, ++n) {
        _lineStream.clear();
        _lineStream.str(*it);

        _lineStream >> k;
        if (_lineStream.fail() || d_pttopo.atomNum(n) != (k - 1)) {
          ostringstream os;
          os << "Error in PERTPOLPARAM block: The atom " << k << " was "
                  << "not found in perturbation topology. Found atom " << d_pttopo.atomNum(n) + 1
                  << " instead.";
          throw InPtTopology::Exception(os.str());
        }

        double alpha[2], dampingLevel[2];
        _lineStream >> l >> nm >> alpha[0] >> dampingLevel[0] >> alpha[1]
                >> dampingLevel[1];

        if (_lineStream.fail()) {
          ostringstream os;
          os << "Bad line in PERTPOLPARAM block (line " << n + 1 << ").";
          throw InPtTopology::Exception(os.str());
        }

        for (unsigned int i = 0; i < 2; ++i) {
          d_pttopo.setPolarisability(n, i, alpha[i]);
          d_pttopo.setDampingLevel(n, i, dampingLevel[i]);
        }
      }
      if (n != numat)
        throw InPtTopology::Exception("Perturbation topology file " + name() +
              " is corrupted. Failed to read all atoms");

      d_pttopo.setHasPolarisationParameters(true);

    }// PERTPOLPARAM
    else {
      d_pttopo.setHasPolarisationParameters(false);
    }
    if (d_blocks.count("MPERTATOM")) {
        std::cerr << "# InPtTopology: WARNING: both PERTATOMPARAM and MPERTATOM block found, ignoring MPERTATOM." << std::endl;
    }
  } else if (d_blocks.count("MPERTATOM")) {
    buffer.clear();
    buffer = d_blocks["MPERTATOM"];
    
    d_pttopo.setMultiPt(true);

    it = buffer.begin() + 1;

    if (buffer.size() < 3)
      throw InPtTopology::Exception("Topology file " + name() +
            " is corrupted. Empty MPERTATOM block!");
    if (buffer[buffer.size() - 1].find("END") != 0)
      throw InPtTopology::Exception("Topology file " + name() +
            " is corrupted. No END in MPERTATOM"
            " block. Got\n"
            + buffer[buffer.size() - 1]);
    _lineStream.clear();
    _lineStream.str(*it);
    _lineStream >> numat >> numpt;
    if (_lineStream.fail()) {
      throw InPtTopology::Exception("Bad line in MPERTATOM block: Could not "
              "read number of atoms or perturbations.");
    }
    d_pttopo.setSize(numat, numpt);
    ++it;
    _lineStream.clear();
    _lineStream.str(*it);
    for (int j = 0; j < d_pttopo.numPt(); ++j) {
      _lineStream >> nm;
      if (_lineStream.fail()) {
        std::ostringstream os;
        os << "Bad line in MPERTATOM block: Could not read perturbation name "
                << j + 1 << "." << endl;
        throw InPtTopology::Exception(os.str());
      }
      d_pttopo.setPertName(j, nm);
    }
    ++it;

    for (n = 0; it < buffer.end() - 1; ++it, ++n) {
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> k >> nm;

      if (n >= numat) {
        std::ostringstream os;
        os << "Bad line in MPERTATOM block: found more atom lines than specified by NJLA"
                << n + 1 << "." << endl;
        throw InPtTopology::Exception(os.str());
      }
      if (_lineStream.fail()) {
        std::ostringstream os;
        os << "Bad line in MPERTATOM block: Could not read atom number or name in line "
                << n + 1 << "." << endl;
        throw InPtTopology::Exception(os.str());
      }

      d_pttopo.setAtomNum(n, k - 1);
      d_pttopo.setAtomName(n, nm);
      for (int j = 0; j < d_pttopo.numPt(); ++j) {
        _lineStream >> iiac[0] >> dq[0];
        if (_lineStream.fail()) {
          std::ostringstream os;
          os << "Bad line in MPERTATOM block: Could not read IAC or charge of state "
                  << j + 1 << " in line " << n + 1 << "." << std::endl;
          throw InPtTopology::Exception(os.str());
        }
        d_pttopo.setIac(n, j, iiac[0] - 1);
        d_pttopo.setCharge(n, j, dq[0]);
      }
      _lineStream >> alphaLJ >> alphaCRF;
      d_pttopo.setAlphaLJ(n, alphaLJ);
      d_pttopo.setAlphaCRF(n, alphaCRF);
    }
    if (n != numat)
      throw InPtTopology::Exception("MPERTATOM block: Perturbation topology file " + name() +
            " is corrupted. Failed to read all atoms");
    if (d_blocks.count("PERTATOMPAIR") || d_blocks.count("PERTBONDSTRETCH") 
    || d_blocks.count("PERTBONDANGLE")|| d_blocks.count("PERTIMPROPERDIH")
    || d_blocks.count("PERTPROPERDIH")|| d_blocks.count("PERTCROSSDIH")) 
        std::cerr << "# InPtTopology: Warning: MPERTATOM block found, ignoring any \n"
                  << "#        additional blocks specifying bonded perturbations." << std::endl;
            
    return; // skip bonded!
  } else {
    throw InPtTopology::Exception("No MPERTATOM or PERTATOMPARAM block in "
            "perturbation topology file\n");
  }

  if (d_blocks.count("PERTATOMPAIR")) {
    int num = _initBlock(buffer, it, "PERTATOMPAIR");
    for (n = 0; it != buffer.end() - 1; ++it, ++n) {
      int i, j, iet[2];
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i >> j >> iet[0] >> iet[1];
      if (_lineStream.fail())
        throw InPtTopology::Exception("Bad line in PERTATOMPAIR block.");

      AtomPairParam apa(i - 1, j - 1, iet[0]), apb(i - 1, j - 1, iet[1]);
      d_pttopo.atompairs(0).insert(apa);
      d_pttopo.atompairs(1).insert(apb);
    }
    if (n != num)
      throw InPtTopology::Exception("Not enough or too many lines in PERTATOMPAIR block.");
  }

  for (unsigned int blockH = 0; blockH <= 1; ++blockH) {
    string block("PERTBONDSTRETCH");
    if (blockH) block += "H";
    if (d_blocks.count(block)) {
      int num = _initBlock(buffer, it, block);
      for (n = 0; it != buffer.end() - 1; ++it, ++n) {
        int i, j, icb[2];
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> i >> j >> icb[0] >> icb[1];
        if (_lineStream.fail())
          throw InPtTopology::Exception("Bad line in " + block + " block.");

        Bond ba(i - 1, j - 1), bb(i - 1, j - 1);
        ba.setType(icb[0] - 1);
        bb.setType(icb[1] - 1);
        d_pttopo.bonds(0).insert(ba);
        d_pttopo.bonds(1).insert(bb);
      }
      if (n != num)
        throw InPtTopology::Exception("Not enough or too many lines in " + block + " block.");
    }
  }

  for (unsigned int blockH = 0; blockH <= 1; ++blockH) {
    string block("PERTBONDANGLE");
    if (blockH) block += "H";
    if (d_blocks.count(block)) {
      int num = _initBlock(buffer, it, block);
      for (n = 0; it != buffer.end() - 1; ++it, ++n) {
        int i, j, k, ica[2];
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> i >> j >> k >> ica[0] >> ica[1];
        if (_lineStream.fail())
          throw InPtTopology::Exception("Bad line in " + block + " block.");

        Angle aa(i - 1, j - 1, k - 1), ab(i - 1, j - 1, k - 1);
        aa.setType(ica[0] - 1);
        ab.setType(ica[1] - 1);
        d_pttopo.angles(0).insert(aa);
        d_pttopo.angles(1).insert(ab);
      }
      if (n != num)
        throw InPtTopology::Exception("Not enough or too many lines in " + block + " block.");
    }
  }

  for (unsigned int blockH = 0; blockH <= 1; ++blockH) {
    string block("PERTIMPROPERDIH");
    if (blockH) block += "H";
    if (d_blocks.count(block)) {
      int num = _initBlock(buffer, it, block);
      for (n = 0; it != buffer.end() - 1; ++it, ++n) {
        int i, j, k, l, ici[2];
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> i >> j >> k >> l >> ici[0] >> ici[1];
        if (_lineStream.fail())
          throw InPtTopology::Exception("Bad line in " + block + " block.");

        Improper ia(i - 1, j - 1, k - 1, l - 1), ib(i - 1, j - 1, k - 1, l - 1);
        ia.setType(ici[0] - 1);
        ib.setType(ici[1] - 1);
        d_pttopo.impropers(0).insert(ia);
        d_pttopo.impropers(1).insert(ib);
      }
      if (n != num)
        throw InPtTopology::Exception("Not enough or too many lines in " + block + " block.");
    }
  }
  for (unsigned int blockH = 0; blockH <= 1; ++blockH) {
    string block("PERTPROPERDIH");
    if (blockH) block += "H";
    if (d_blocks.count(block)) {
      int num = _initBlock(buffer, it, block);
      for (n = 0; it != buffer.end() - 1; ++it, ++n) {
        int i, j, k, l, icd[2];
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> i >> j >> k >> l >> icd[0] >> icd[1];
        if (_lineStream.fail())
          throw InPtTopology::Exception("Bad line in " + block + " block.");

        Dihedral da(i - 1, j - 1, k - 1, l - 1), db(i - 1, j - 1, k - 1, l - 1);
        da.setType(icd[0] - 1);
        db.setType(icd[1] - 1);
        d_pttopo.dihedrals(0).insert(da);
        d_pttopo.dihedrals(1).insert(db);
      }
      if (n != num)
        throw InPtTopology::Exception("Not enough or too many lines in " + block + " block.");
    }
  }
  for (unsigned int blockH = 0; blockH <= 1; ++blockH) {
    string block("PERTCROSSDIH");
    if (blockH) block += "H";
    if (d_blocks.count(block)) {
      int num = _initBlock(buffer, it, block);
      for (n = 0; it != buffer.end() - 1; ++it, ++n) {
        int a, b, c, d, e, f, g, h, icd[2];
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> a >> b >> c >> d >> e >> f >> g >> h >> icd[0] >> icd[1];
        if (_lineStream.fail())
          throw InPtTopology::Exception("Bad line in " + block + " block.");

        CrossDihedral da(a - 1, b - 1, c - 1, d - 1, e - 1, f - 1, g - 1, h - 1), db(a - 1, b - 1, c - 1, d - 1, e - 1, f - 1, g - 1, h - 1);
        da.setType(icd[0] - 1);
        db.setType(icd[1] - 1);
        d_pttopo.crossdihedrals(0).insert(da);
        d_pttopo.crossdihedrals(1).insert(db);
      }
      if (n != num)
        throw InPtTopology::Exception("Not enough or too many lines in " + block + " block.");
    }
  }
}





