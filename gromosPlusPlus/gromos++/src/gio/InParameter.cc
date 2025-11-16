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

// gio_InParameter.cc
#include "InParameter.h"

#include <cassert>
#include <cmath>
#include <map>
#include <map>
#include <vector>
#include <string>

#include "Ginstream.h"
#include "../gcore/MassType.h"
#include "../gcore/VirtualAtomType.h"
#include "../gcore/BondType.h"
#include "../gcore/AngleType.h"
#include "../gcore/DihedralType.h"
#include "../gcore/ImproperType.h"
#include "../gcore/LJType.h"
#include "../gcore/LJExceptionType.h"
#include "../gcore/AtomPair.h"
#include "../gcore/GromosForceField.h"
#include "../gmath/Physics.h"
#include "../args/Arguments.h"


using namespace std;
using namespace gcore;

// Implementation class
class gio::InParameter_i: public gio::Ginstream
{
  friend class gio::InParameter;
  gcore::GromosForceField d_gff;
  std::map<std::string, std::vector<std::string> > d_blocks;
  /**
   * The init function reads in the whole file into the map of blocks.
   */
  void init();
  /**
   * _initBlock is a function that initializes the reading of
   * a block. It checks whether the block is read in, and returns
   * the number of lines that are to be read in from it.
   */
  int _initBlock(std::vector<std::string> &buffer,
		 std::vector<std::string>::const_iterator &it,
		 const string blockname);
  /**
   * parseForceField takes all relevant blocks and stores the information in
   * d_gff
   */
  void parseForceField();
  
  InParameter_i (std::string &s): d_gff(), d_blocks()
  {
    this->open(s);
    this->init();
    this->parseForceField();
  }
};

// Constructors

gio::InParameter::InParameter(std::string name){
  d_this = new InParameter_i(name);
}

gio::InParameter::~InParameter(){
  delete d_this;
}

const std::string gio::InParameter::title()const{
  return d_this->title();
}

const gcore::GromosForceField &gio::InParameter::forceField()const{
  return d_this->d_gff;
}

int gio::InParameter_i::_initBlock(std::vector<std::string> &buffer,
				  std::vector<std::string>::const_iterator &it,
				   const string blockname)
{
  buffer.clear();
  buffer=d_blocks[blockname];
  if(buffer.size() < 2)
    throw InParameter::Exception("Parameter file "+name()+
				" is corrupted. No "+blockname+
				" block!");
  if(buffer[buffer.size()-1].find("END")!=0)
      throw InParameter::Exception("Topology file " + name() +
				       " is corrupted. No END in "+blockname+
				       " block. Got\n"
				       + buffer[buffer.size()-1]);

  it=buffer.begin()+1;
  return buffer.size()-2;
}

void gio::InParameter_i::init(){

  if(!stream())
    throw InParameter::Exception("Could not open parameter file "+name());

  // First read the whole file into the map
  std::vector<std::string> buffer;
  
  while(!stream().eof()){
    getblock(buffer);
    if(buffer.size()){
      d_blocks[buffer[0]] = buffer;
      buffer.clear();
    }
  }
}


void gio::InParameter_i::parseForceField()
{
  // generic variables
  double d[4];
  int num, n, i[5];
  string s;
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  { // FORCEFIELD block
    buffer.clear();
    buffer=d_blocks["FORCEFIELD"];
    if(buffer.size()>0){
      if(buffer.size() != 3) 
	throw InParameter::Exception("Parameter file " + name() +
				     " is corrupted. FORCEFIELD block should have only "
				     "one line");
      if(buffer[buffer.size()-1].find("END")!=0)
	throw InParameter::Exception("Parameter file " + name() +
				     " is corrupted. No END in FORCEFIELD"
				     " block. Got\n"
				     + buffer[buffer.size()-1]);

      d_gff.setForceField(buffer[1]);
    }
  }
  

  if (args::Arguments::inG96 == true) {
    { // MASSATOMTYPECODE block
      num = _initBlock(buffer, it, "MASSATOMTYPECODE");
      for (n = 0; n < num; ++it, ++n) {
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> i[0] >> d[0] >> s;
        if (_lineStream.fail())
          throw InParameter::Exception("Bad line in MASSATOMTYPECODE block:\n"
                + *it);
        d_gff.addMassType(MassType(--i[0], d[0]));
      }
    } // MASSATOMTYPECODE block
  } else {
    { // MASSATOMTYPECODE block
      num = _initBlock(buffer, it, "MASSATOMTYPECODE");
      int num_mass_type_codes;
      int largest_mass_type_code;
      // Read in NRMATY and NMATY
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> num_mass_type_codes >> largest_mass_type_code;
      if (_lineStream.fail())
        throw InParameter::Exception("Bad line in MASSATOMTYPECODE block:\n"
              + *it);
      if ((num - 1) != num_mass_type_codes)
        throw InParameter::Exception("Not enough or too many lines in "
              "MASSATOMTYPECODE block:");

      ++it;
      for (n = 0; n < num-1; ++it, ++n) {
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> i[0] >> d[0] >> s;
        if (_lineStream.fail())
          throw InParameter::Exception("Bad line in MASSATOMTYPECODE block:\n"
                + *it);
        
        if (i[0] > largest_mass_type_code)
        throw InParameter::Exception(
              "MASSATOMTYPECODE block: Mass type code larger than maximum (NMATY)");
        
        d_gff.addMassType(MassType(--i[0], d[0]));
      }
    } // MASSATOMTYPECODE block    
  }
  { // VIRTUALATOMTYPECODE block
    buffer=d_blocks["VIRTUALATOMTYPECODE"];
    if(buffer.size()>0){
      num = _initBlock(buffer, it, "VIRTUALATOMTYPECODE");
      int num_va_type_codes;
      int largest_va_type_code;
      // Read in NRBTY and NBTY
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> num_va_type_codes >> largest_va_type_code;
      if(_lineStream.fail())
        throw InParameter::Exception("Bad line in VIRTUALATOMTYPECODE block:\n"
                                             +*it);
      if ((num-1) != num_va_type_codes)
        throw InParameter::Exception("Not enough or too many lines in "
                "VIRTUALATOMTYPECODE block:");

      ++it;

      for(n=0; n<num-1; ++it, ++n){
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> i[0] >> d[0] >> d[1];
        if (_lineStream.fail())
          throw InParameter::Exception("Bad line in VIRTUALATOMTYPECODE block:\n"
                + *it);

        if (i[0] > largest_va_type_code)
          throw InParameter::Exception(
                "VIRTUALATOMTYPECODE block: Virtual atom type code larger than maximum (NVATY)");

        d_gff.addVirtualAtomType(VirtualAtomType(i[0], d[0], d[1]));
      }
    }
  } // VIRTUALATOMTYPECODE block
  if(args::Arguments::inG96==true){ // BONDTYPECODE block
    num = _initBlock(buffer, it, "BONDTYPECODE");
    for(n=0; n<num; ++it, ++n){
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i[0] >> d[0] >> d[1];
      if(_lineStream.fail())
	throw InParameter::Exception("Bad line in BONDTYPECODE block:\n"
				     +*it);
      if(i[0]!=n+1)
	throw InParameter::Exception(
	     "BondTypes in BONDTYPECODE block are not sequential");
     
      // Convert quartic force constant into quadratic one (Vol. 5, p. 10) 
      d[2] = 2.0 * d[1] * d[1] *d[0];
      d_gff.addBondType(BondType(--i[0], d[0], d[1]));
    }
  } // BONDTYPECODE block (g96 only)
  else{ // BONDSTRETCHTYPECODE block
    num = _initBlock(buffer, it, "BONDSTRETCHTYPECODE");
    int num_bond_type_codes;
    int largest_bond_type_code;
    // Read in NRBTY and NBTY
    _lineStream.clear();
    _lineStream.str(*it);
    _lineStream >> num_bond_type_codes >> largest_bond_type_code;
    if(_lineStream.fail())
      throw InParameter::Exception("Bad line in BONDSTRETCHTYPECODE block:\n"
                                             +*it);
    if ((num-1) != num_bond_type_codes) 
      throw InParameter::Exception("Not enough or too many lines in "
              "BONDSTRETCHTYPECODE block:");
    
    ++it;

    for(n=0; n<num-1; ++it, ++n){
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i[0] >> d[0] >> d[1] >> d[2];
      if (_lineStream.fail())
        throw InParameter::Exception("Bad line in BONDSTRETCHTYPECODE block:\n"
              + *it);
      
      if (i[0] > largest_bond_type_code)
        throw InParameter::Exception(
              "BONDSTRETCHTYPECODE block: Bond type code larger than maximum (NBTY)");

      d_gff.addBondType(BondType(--i[0], d[0], d[1], d[2]));
    }
  } // BONDSTRETCHTYPECODE block
  if(args::Arguments::inG96==true){ // BONDANGLETYPECOD block
    num = _initBlock(buffer, it, "BONDANGLETYPECOD");
    for(n=0; n<num; ++it, ++n){
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i[0] >> d[0] >> d[1];
      if(_lineStream.fail())
	throw InParameter::Exception("Bad line in BONDANGLETYPECOD block:\n"
				     +*it);
      if(i[0]!=n+1)
	throw InParameter::Exception(
	     "AngleTypes in BONDANGLETYPECOD block are not sequential");
  
      try {
        d_gff.addAngleType(AngleType(i[0]-1, d[0], d[1]));
      } catch (gromos::Exception & exp) {
        if (!args::Arguments::outG96) {
          std::cerr << exp.what() << std::endl
                  << "Setting harmonic force constant to -1.0." << std::endl;
        }
        d_gff.addAngleType(AngleType(i[0]-1, d[0], -1.0, d[1]));
      }
    }
  } //BONDANGLETYPECODE
  else { // BONDANGLEBENDTYPECODE block
    num = _initBlock(buffer, it, "BONDANGLEBENDTYPECODE");
    int num_angle_type_codes;
    int largest_angle_type_code;
    // Read in NRTTY and NTTY
    _lineStream.clear();
    _lineStream.str(*it);
    _lineStream >> num_angle_type_codes >> largest_angle_type_code;
    if(_lineStream.fail())
      throw InParameter::Exception("Bad line in BONDANGLEBENDTYPECODE block:\n"
                                             +*it);
    if ((num-1) != num_angle_type_codes) 
      throw InParameter::Exception("Not enough or too many lines in "
              "BONDANGLEBENDTYPECODE block:");
    
    ++it; 

    for(n=0; n<num-1; ++it, ++n) {
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i[0] >> d[0] >> d[1] >> d[2];
      if (_lineStream.fail())
        throw InParameter::Exception("Bad line in BONDANGLEBENDTYPECODE block:\n"
              + *it);
      if (i[0] > largest_angle_type_code)
        throw InParameter::Exception(
              "BONDSTRETCHTYPECODE block: Angle type code larger than maximum (NPTY)");

      d_gff.addAngleType(AngleType(--i[0], d[0], d[1], d[2]));
    }
  } //BONDANGLEBENDTYPECODE
  if(args::Arguments::inG96==true){ // IMPDIHEDRALTYPEC block
    num = _initBlock(buffer, it, "IMPDIHEDRALTYPEC");
    for(n=0; n<num; ++it, ++n){
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i[0] >> d[0] >> d[1];
      if(_lineStream.fail())
	throw InParameter::Exception("Bad line in IMPDIHEDRALTYPEC block:\n"
				     +*it);
      if(i[0]!=n+1)
	throw InParameter::Exception(
	     "ImproperTypes in IMPDIHEDRALTYPEC block are not sequential");
      d_gff.addImproperType(ImproperType(--i[0], d[0],d[1]));
    }
  } // IMPDIHEDRALTYPEC
  else { // IMPDIHEDRALTYPECODE block
    num = _initBlock(buffer, it, "IMPDIHEDRALTYPECODE");
    // Read in NRQTY and NQTY
    int num_impdihedral_type_codes;
    int largest_impdihedral_type_code;
    _lineStream.clear();
    _lineStream.str(*it);
    _lineStream >> num_impdihedral_type_codes >> largest_impdihedral_type_code;
    if(_lineStream.fail())
      throw InParameter::Exception("Bad line in IMPDIHEDRALTYPECODE block:\n"
                                             +*it);
    if ((num-1) != num_impdihedral_type_codes) 
      throw InParameter::Exception("Not enough or too many lines in "
              "IMPDIHEDRALTYPECODE block:");
    ++it;

    for(n=0; n<num-1; ++it, ++n){
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i[0] >> d[0] >> d[1];
      if (_lineStream.fail())
        throw InParameter::Exception("Bad line in IMPDIHEDRALTYPECODE block:\n"
              + *it);
      if (i[0] > largest_impdihedral_type_code)
        throw InParameter::Exception(
              "BONDSTRETCHTYPECODE block: Improper dihedral type code larger "
              "than maximum (NQTY)");
      d_gff.addImproperType(ImproperType(--i[0], d[0], d[1]));
    }
  } // IMPDIHEDRALTYPECODE
  if(args::Arguments::inG96==true){ // DIHEDRALTYPECODE block
    num = _initBlock(buffer, it, "DIHEDRALTYPECODE");
    for(n=0; n<num; ++it, ++n){
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i[0] >> d[0] >> d[1] >> i[1];
      if(_lineStream.fail())
	throw InParameter::Exception("Bad line in DIHEDRALTYPECODE block:\n"
				     +*it);
      if(i[0]!=n+1)
	throw InParameter::Exception(
	     "DihedralTypes in DIHEDRALTYPECODE block are not sequential");
      // Convert phase into phase-shift angle(given in degrees)
      d[2] = acos(d[1]) * gmath::physConst.get_radian2degree();
      d_gff.addDihedralType(DihedralType(--i[0], d[0], d[1], d[2], i[1]));
    }
  } // DIHEDRALTYPECODE
  else { // TORSDIHEDRALTYPECODE block
    num = _initBlock(buffer, it, "TORSDIHEDRALTYPECODE");
    // Read in NRPTY and NPTY
    int num_dihedral_type_codes;
    int largest_dihedral_type_code;
    _lineStream.clear();
    _lineStream.str(*it);
    _lineStream >> num_dihedral_type_codes >> largest_dihedral_type_code;
    if(_lineStream.fail())
      throw InParameter::Exception("Bad line in TORSDIHEDRALTYPECODE block:\n"
                                             +*it);
    if ((num-1) != num_dihedral_type_codes) 
      throw InParameter::Exception("Not enough or too many lines in "
              "TORSDIHEDRALTYPECODE block:");
    ++it;

    for(n=0; n<num-1; ++it, ++n){
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i[0] >> d[0] >> d[2] >> i[1];
      if (_lineStream.fail())
        throw InParameter::Exception("Bad line in TORSDIHEDRALTYPECODE block:\n"
              + *it);
      if (i[0] > largest_dihedral_type_code)
        throw InParameter::Exception(
              "TORSDIHEDRALTYPECODE block: dihedral type code larger than "
              "maximum (NPTY)");

      // Convert phase-shift angle(given in degrees) into phase
      d[1] = cos(d[2]*gmath::physConst.get_degree2radian());
      d_gff.addDihedralType(DihedralType(--i[0], d[0], d[1], d[2], i[1]));
    }
  } // TORSDIHEDRALTYPECODE
  { // SINGLEATOMLJPAIR
    num = _initBlock(buffer, it, "SINGLEATOMLJPAIR");
    _lineStream.clear();
    _lineStream.str(*it);
    _lineStream >> num;
    std::vector<double> sc6(num), sc12[3], scs6(num), scs12(num);
    sc12[0].resize(num);
    sc12[1].resize(num);
    sc12[2].resize(num);
    
    std::vector<std::vector<int> >    pl(num);
    for(int j=0; j<num; ++j){
      pl[j].resize(num);
    }
    std::string ljblock;
    gio::concatenate(it+1, buffer.end()-1, ljblock);
    _lineStream.clear();
    _lineStream.str(ljblock);
    for(n=0; n<num; n++){
      _lineStream >> i[0]>> s>> sc6[n]>> sc12[0][n]>> sc12[1][n]>> sc12[2][n];
      _lineStream >> scs6[n] >> scs12[n];
      if(_lineStream.fail()){
	ostringstream os;
	os << "Bad line in SINGLEATOMLJPAIR block, IAC: " << n+1 << "\n"
	   << "Trying to read parameters";
	throw InParameter::Exception(os.str());
      }
      if(i[0]!=n+1)
	throw InParameter::Exception(
	     "AtomTypes in SINGLEATOMLJPAIR block are not sequential");

      for(int k=0; k<num; k++)
	_lineStream >> pl[n][k];
      if(_lineStream.fail()){
	ostringstream os;
	os << "Bad line in SINGLEATOMLJPAIR block, IAC: " << n+1 << "\n"
	   << "Trying to read " << num << " elements of the interaction "
	   << "matrix";
	throw InParameter::Exception(os.str());
      }
      d_gff.addAtomTypeName(s);
      for(int k=0; k<=n; k++){
	d[1] = sc6[n]              * sc6[k];
	d[0] = sc12[pl[n][k]-1][n] * sc12[pl[k][n]-1][k];
	d[3] = scs6[n]             * scs6[k];
	d[2] = scs12[n]            * scs12[k];
	d_gff.setLJType(AtomPair(n,k), LJType(d[0], d[1], d[2], d[3]));
      }
    }
  } // SINGLEATOMLJPAIR
  { // MIXEDATOMLJPAIR block
    // this one does not have to be there
    if(d_blocks.count("MIXEDATOMLJPAIR")){
      
      num = _initBlock(buffer, it, "MIXEDATOMLJPAIR");
      for(n =0; n<num; ++it, ++n){
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> i[0] >> i[1] >> d[1] >> d[0] >> d[3] >> d[2];
	if(_lineStream.fail())
	  throw InParameter::Exception("Bad line in MIXEDATOMLJPAIR block:\n"
				       +*it);
	d_gff.setLJType(AtomPair(--i[0],--i[1]),LJType(d[0],d[1],d[2],d[3]));
      }
    } // MIXEDATOMLJPAIR
  }
  { // SPECATOMLJPAIR block
    num = _initBlock(buffer, it, "SPECATOMLJPAIR");
    for (n = 0; n < num; ++it, ++n) {
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i[0] >> d[0] >> d[1];
      if (_lineStream.fail())
        throw InParameter::Exception("Bad line in SPECATOMLJPAIR block:\n" + *it);
      d_gff.addLJExceptionType(LJExceptionType(--i[0], d[0], d[1]));
    }
  } // end SPECATOMLJPAIR block
  
}




