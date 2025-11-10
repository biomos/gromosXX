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
#include "InPDB.h"

#include <algorithm>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cassert>
#include <iostream>
#include <iterator>

#include "../gromos/Exception.h"
#include "../gmath/Vec.h"

// the pdb file format definitions (positions of data information)
// for more information see:
// http://www.wwpdb.org/documentation/format32/sect9.html
#define TYPE_POS 0
#define TYPE_LEN 6
#define TYPE substr(TYPE_POS,TYPE_LEN)
#define SERIAL_POS 6
#define SERIAL_LEN 5
#define SERIAL substr(SERIAL_POS,SERIAL_LEN)
#define ATOMNAME_POS 12
#define ATOMNAME_LEN 4
#define ATOMNAME substr(ATOMNAME_POS, ATOMNAME_LEN)
// In the official PDB this is a separate "column", actually...
//#define ALTLOC_POS 16
//#define ALTLOC_LEN 1
//#define ALTLOC substr(ALTLOC_POS,ALTLOC_LEN)
#define RESNAME_POS 16
#define RESNAME_LEN 4 // this is different from the official PDB format but
                      // needed for GROMOS
#define RESNAME substr(RESNAME_POS, RESNAME_LEN)
#define CHAINID_POS 21
#define CHAINID_LEN 1
#define CHAINID substr(CHAINID_POS,CHAINID_LEN)
#define RESNUM_POS 22
#define RESNUM_LEN 4
#define RESNUM substr(RESNUM_POS, RESNUM_LEN)
#define ICODE_POS 26
#define ICODE_LEN 1
#define ICODE substr(ICODE_POS, ICODE_LEN)
#define COORDX_POS 30
#define COORDX_LEN 8
#define COORDX substr(COORDX_POS, COORDX_LEN)
#define COORDY_POS 38
#define COORDY_LEN 8
#define COORDY substr(COORDY_POS, COORDY_LEN)
#define COORDZ_POS 46
#define COORDZ_LEN 8
#define COORDZ substr(COORDZ_POS, COORDZ_LEN)
#define OCCUPANCY_POS 54
#define OCCUPANCY_LEN 6
#define OCCUPANCY substr(OCCUPANCY_POS,OCCUPANCY_LEN)
#define BFACTOR_POS 60
#define BFACTOR_LEN 6
#define BFACTOR substr(BFACTOR_POS,BFACTOR_LEN)
#define ELEMENT_POS 76
#define ELEMENT_LEN 2
#define ELEMENT substr(ELEMENT_POS,ELEMENT_LEN)
#define CHARGE_POS 78
#define CHARGE_LEN 2
#define CHARGE substr(CHARGE_POS, CHARGE_POS)

std::string removeChar(std::string s, char c = ' ');

using namespace std;

namespace gio {

  /**
   * Class InPDB_ is the implementations class of the class PDB. It holds
   * all the member data so the PDB class needs nothing more than a pointer
   * to this class for its data.
   */
  class InPDB_i {
  public:

    /**
     * The file name of the pdb file which will be read.
     */
    string filemane;
    /**
     * A switch to read or not to read the ATOM atoms of the pdb file.
     */
    bool readATOM;
    /**
     * A switch to read or not to read the HETATM atoms of the pdb file.
     */
    bool readHETATM;
    /**
     * A vector to store the atom types of the pdb atoms.
     */
    vector<string> types;
    /**
     * A vector to store the serial number of the pdb atoms.
     */
    vector<int> serials;
    /**
     * A vector to store the atom names of the pdb atoms.
     */
    vector<string> atoms;
    /**
     * A vector to store the residue names of the pdb atoms.
     */
    vector<string> resNames;
    /**
     * A vector to store the chain IDs of the pdb atoms.
     */
    vector<char> chainIDs;
    /**
     * A vector to store the chain IDs of the pdb atoms.
     */
    vector<int> resNums;
    /**
     * A vector to store the insertion codes.
     */
    vector<char> iCodes;
    /**
     * A vector to store the x-coordinates of the pdb atoms.
     */
    vector<double> X;
    /**
     * A vector to store the y-coordinates of the pdb atoms.
     */
    vector<double> Y;
    /**
     * A vector to store the z-coordinates of the pdb atoms.
     */
    vector<double> Z;
    /**
     * A vector to store the occupancies of the pdb atoms.
     */
    vector<double> occupancies;
    /**
     * A vector to store the temp factors of the pdb atoms.
     */
    vector<double> tempFactors;
    /**
     * A vector to store the element symbols of the pdb atoms.
     */
    vector<string> elements;
    /**
     * A vector to store the atom charges of the pdb atoms.
     */
    vector<string> charges;
    /**
     * A vector to store the charges of the pqr atoms.
     */
    vector<double> pqr_charges;
    /**
     * A vector to store the radii of the pqr atoms.
     */
    vector<double> pqr_radii;
    
    bool b_type;
    bool b_serial;
    bool b_atomname;
    bool b_resname;
    bool b_chainid;
    bool b_resseq;
    bool b_icode;
    bool b_x;
    bool b_y;
    bool b_z;
    bool b_occupancy;
    bool b_tempfactor;
    bool b_element;
    bool b_charge;
    
    /**
     * A vector to store the residue sequences as it is in the pdb file.
     */
    vector<string> resSeq;

  };

  // Constructor

  InPDB::InPDB(const std::string &filename, bool readATOM, bool readHETATM) {
    d_this = new InPDB_i;
    d_this->filemane = filename;
    d_this->readATOM = readATOM;
    d_this->readHETATM = readHETATM;
    // the following booleans indicate if the corresponding variable is
    // specified in the pdb file or not
    d_this->b_serial = false;
    d_this->b_atomname = false;
    d_this->b_resname = false;
    d_this->b_chainid = false;
    d_this->b_resseq = false;
    d_this->b_icode = false;
    d_this->b_x = false;
    d_this->b_y = false;
    d_this->b_z = false;
    d_this->b_occupancy = false;
    d_this->b_tempfactor = false;
    d_this->b_element = false;
    d_this->b_charge = false;
  }

  InPDB::~InPDB() {
    if (d_this)delete d_this;
  }

  void InPDB::select(const std::string &thing) {
    if (thing == "ATOM") {
      d_this->readATOM = true;
      d_this->readHETATM = false;
    } else if (thing == "HETATM") {
      d_this->readATOM = false;
      d_this->readHETATM = true;
    } else if (thing == "ALL") {
      d_this->readATOM = true;
      d_this->readHETATM = true;
    } else {
      stringstream msg;
      msg << "InPDB::select does not know the argument " << thing;
      throw Exception(msg.str());
    }
  }

  void InPDB::read() {

    // open the PDB file
    ifstream fin(d_this->filemane.c_str());

    // check if the PDB file could be opened
    if (!fin.is_open()) {
      stringstream msg;
      msg << "Could not open the PDB file " << d_this->filemane;
      throw InPDB::Exception(msg.str());
    }

    // read and save the contents of the PDB file
    string line;
    stringstream ssline;
    int res = -1;
    while (!fin.eof()) {
      getline(fin, line);
      string type = "";
      if (line.size() > TYPE_POS) {
        ssline << line.TYPE << endl;
        ssline >> type;
      }
      ssline.clear();
      ssline.str("");
      if ((type == "ATOM" && d_this->readATOM) ||
	  (type == "HETATM" && d_this->readHETATM)) {
        if (removeChar(line.SERIAL) != "") {
          d_this->b_serial = true;
          ssline << line.SERIAL << endl;
        } else if (d_this->b_serial) {
          stringstream msg;
          msg << "bad line in pdb file (" << d_this->filemane << "):\n" << line << endl;
          throw InPDB::Exception(msg.str());
        }
        if (removeChar(line.ATOMNAME) != "") {
          d_this->b_atomname = true;
          ssline << line.ATOMNAME << endl;
        } else if (d_this->b_atomname) {
          stringstream msg;
          msg << "bad line in pdb file (" << d_this->filemane << "):\n" << line << endl;
          throw InPDB::Exception(msg.str());
        }
        if (removeChar(line.RESNAME) != "") {
          d_this->b_resname = true;
          ssline << line.RESNAME << endl;
        } else if (d_this->b_resname) {
          stringstream msg;
          msg << "bad line in pdb file (" << d_this->filemane << "):\n" << line << endl;
          throw InPDB::Exception(msg.str());
        }
        if (removeChar(line.CHAINID) != "") {
          d_this->b_chainid = true;
          ssline << line.CHAINID << endl;
        } else if (d_this->b_chainid) {
          stringstream msg;
          msg << "bad line in pdb file (" << d_this->filemane << "):\n" << line << endl;
          throw InPDB::Exception(msg.str());
        }
        if (removeChar(line.RESNUM) != "") {
          d_this->b_resseq = true;
          ssline << line.RESNUM << endl;
        } else if (d_this->b_resseq) {
          stringstream msg;
          msg << "bad line in pdb file (" << d_this->filemane << "):\n" << line << endl;
          throw InPDB::Exception(msg.str());
        }
        if (removeChar(line.ICODE) != "") {
          d_this->b_icode = true;
          ssline << line.ICODE << endl;
        } else if (d_this->b_icode) {
          stringstream msg;
          msg << "bad line in pdb file (" << d_this->filemane << "):\n" << line << endl;
          throw InPDB::Exception(msg.str());
        }
        if (removeChar(line.COORDX) != "") {
          d_this->b_x = true;
          ssline << line.COORDX << endl;
        } else if (d_this->b_x) {
          stringstream msg;
          msg << "bad line in pdb file (" << d_this->filemane << "):\n" << line << endl;
          throw InPDB::Exception(msg.str());
        }
        if (removeChar(line.COORDY) != "") {
          d_this->b_y = true;
          ssline << line.COORDY << endl;
        } else if (d_this->b_y) {
          stringstream msg;
          msg << "bad line in pdb file (" << d_this->filemane << "):\n" << line << endl;
          throw InPDB::Exception(msg.str());
        }
        if (removeChar(line.COORDZ) != "") {
          d_this->b_z = true;
          ssline << line.COORDZ << endl;
        } else if (d_this->b_z) {
          stringstream msg;
          msg << "bad line in pdb file (" << d_this->filemane << "):\n" << line << endl;
          throw InPDB::Exception(msg.str());
        }
        if (removeChar(line.OCCUPANCY) != "") {
          d_this->b_occupancy = true;
          ssline << line.OCCUPANCY << endl;
        } else if (d_this->b_occupancy) {
          stringstream msg;
          msg << "bad line in pdb file (" << d_this->filemane << "):\n" << line << endl;
          throw InPDB::Exception(msg.str());
        }
        if (removeChar(line.BFACTOR) != "") {
          d_this->b_tempfactor = true;
          ssline << line.BFACTOR << endl;
        } else if (d_this->b_tempfactor) {
          stringstream msg;
          msg << "bad line in pdb file (" << d_this->filemane << "):\n" << line << endl;
          throw InPDB::Exception(msg.str());
        }

	if (line.length() > 75) {
	  if (removeChar(line.ELEMENT) != "") {
	    d_this->b_element = true;
	    ssline << line.ELEMENT << endl;
	  } else if (d_this->b_element) {
	    stringstream msg;
	    msg << "bad line in pdb file (" << d_this->filemane << "):\n" << line << endl;
	    throw InPDB::Exception(msg.str());
	  }
	}
	if (line.length() > 77) {
	  if (removeChar(line.CHARGE) != "") {
	    d_this->b_charge = true;
	    ssline << line.CHARGE << endl;
	  } else if (d_this->b_charge) {
	    stringstream msg;
	    msg << "bad line in pdb file (" << d_this->filemane << "):\n" << line << endl;
	    throw InPDB::Exception(msg.str());
	  }
	}

        int serial;
        string atom;
        string resName;
        char chainID;
        int seqNo;
        char iCode;
        double x;
        double y;
        double z;
        double occupancy;
        double tempFactor;
        string element;
        string charge;
        if(d_this->b_serial) {
          ssline >> serial;
        }
        if(d_this->b_atomname) {
          ssline >> atom;
        }
        if(d_this->b_resname) {
          ssline >> resName;
        }
        if(d_this->b_chainid) {
          ssline >> chainID;
        }
        if(d_this->b_resseq) {
          ssline >> seqNo;
        }
        if(d_this->b_icode) {
          ssline >> iCode;
        }
        if(d_this->b_x) {
          ssline >> x;
        }
        if(d_this->b_y) {
          ssline >> y;
        }
        if(d_this->b_z) {
          ssline >> z;
        }
        if(d_this->b_occupancy) {
          ssline >> occupancy;
        }
        if(d_this->b_tempfactor) {
          ssline >> tempFactor;
        }
        if(d_this->b_element) {
          ssline >> element;
        }
        if(d_this->b_charge) {
          ssline >> charge;
        }
        // error message in case of the conversion failed...
        if (ssline.bad() || ssline.fail()) {
          stringstream msg;
          msg << "bad line in PDB file:\n" << line;
          cerr << "ssline.bad() = " << ssline.bad() << ", ssline.fail() = " << ssline.fail() << endl;
          throw gromos::Exception("InPDB.cc", msg.str());
        }
        
        // memorize the read variables
        d_this->types.push_back(type);
        if(d_this->b_serial) {
          d_this->serials.push_back(serial);
        }
        if(d_this->b_atomname) {
          d_this->atoms.push_back(atom);
        }
        if(d_this->b_resname) {
          d_this->resNames.push_back(resName);
        }
        if(d_this->b_chainid) {
          d_this->chainIDs.push_back(chainID);
        }
        if(d_this->b_resseq) {
          d_this->resNums.push_back(seqNo);
        }
        if(d_this->b_icode) {
          d_this->iCodes.push_back(iCode);
        }
        if(d_this->b_x) {
          d_this->X.push_back(x);
        }
        if(d_this->b_y) {
          d_this->Y.push_back(y);
        }
        if(d_this->b_z) {
          d_this->Z.push_back(z);
        }
        if(d_this->b_occupancy) {
          d_this->occupancies.push_back(occupancy);
        }
        if(d_this->b_tempfactor) {
          d_this->tempFactors.push_back(tempFactor);
        }
        if(d_this->b_element) {
          d_this->elements.push_back(element);
        }
        if(d_this->b_charge) {
          d_this->charges.push_back(charge);
        }
        
        // add the residue to the sequence (in case of a new residue)
        if (res != seqNo) {
          res = seqNo;
          d_this->resSeq.push_back(resName);
        }
        ssline.clear();
        ssline.str("");
      }
    } // while
    // close the PDB file
    fin.close();

  } // end of read

  std::vector<std::string> split(std::string const &input) {
    std::istringstream buffer(input);
    std::vector<std::string> ret;
    
    std::copy(std::istream_iterator<std::string>(buffer), 
              std::istream_iterator<std::string>(),
              std::back_inserter(ret));
    return ret;
  }
  
  void InPDB::readPQR() {

    // open the PQR file
    ifstream fin(d_this->filemane.c_str());

    // check if the PQR file could be opened
    if (!fin.is_open()) {
      stringstream msg;
      msg << "Could not open the PQR file " << d_this->filemane;
      throw InPDB::Exception(msg.str());
    }

    // read and save the contents of the PQR file
    string line;
    istringstream issline(line);
    stringstream ssline;
    int res = -1;
    unsigned int elementCount = 0; // for counting the number of elements per line
    std::vector<std::string> lineSplit;
    while (!fin.eof()) {
      getline(fin, line);
      string type = "";
      if (line.size() > TYPE_POS) {
        ssline << line.TYPE << endl;
        ssline >> type;
      }
      ssline.clear();
      ssline.str("");
      if ((type == "ATOM" && d_this->readATOM) ||
	  (type == "HETATM" && d_this->readHETATM)) {


	std::istringstream iss(line);
	for(std::string s; iss >> s; ) {
	  ssline << s << endl;
	  elementCount++;
	}
	
	if (elementCount < 10) {
	  stringstream msg;
	  msg << "At least one ATOM or HETATOM line has less than 10 elements. Don't know what to do with it!" << endl;
	  throw InPDB::Exception(msg.str());
	} else if (elementCount > 11) {
	  stringstream msg;
	  msg << "At least one ATOM or HETATOM line has more than 11 elements. Don't know what to do with it!" << endl;
	  throw InPDB::Exception(msg.str());
	}

	// read the stringstream into variables
	string fieldName;
	int serial;
        string atom;
        string resName;
        char chainID;
        int seqNo;
        double x;
        double y;
        double z;
        double charge;
        double radius;
	if (elementCount == 10) {
	  ssline >> fieldName;
	  ssline >> serial;
	  ssline >> atom;
	  ssline >> resName;
	  ssline >> seqNo;
	  ssline >> x;
	  ssline >> y;
	  ssline >> z;
	  ssline >> charge;
	  ssline >> radius;
	} else if (elementCount == 11) {
	  ssline >> fieldName;
	  ssline >> serial;
	  ssline >> atom;
	  ssline >> resName;
	  ssline >> chainID;
	  ssline >> seqNo;
	  ssline >> x;
	  ssline >> y;
	  ssline >> z;
	  ssline >> charge;
	  ssline >> radius;
	}
	
        // error message in case of the conversion failed...
        if (ssline.bad() || ssline.fail()) {
          stringstream msg;
          msg << "bad line in PDB file:\n" << line;
          cerr << "ssline.bad() = " << ssline.bad() << ", ssline.fail() = " << ssline.fail() << endl;
          throw gromos::Exception("InPDB.cc", msg.str());
        }
        
        // memorize the read variables
	d_this->serials.push_back(serial);
	d_this->atoms.push_back(atom);
	d_this->resNames.push_back(resName);
	if (elementCount == 9) {
	  d_this->chainIDs.push_back(chainID);
	}
	d_this->resNums.push_back(seqNo);
	d_this->X.push_back(x);
	d_this->Y.push_back(y);
	d_this->Z.push_back(z);
	d_this->pqr_charges.push_back(charge);
	d_this->pqr_radii.push_back(radius);
	  
        
        // add the residue to the sequence (in case of a new residue)
        if (res != seqNo) {
          res = seqNo;
          d_this->resSeq.push_back(resName);
        }
        ssline.clear();
        ssline.str("");
	elementCount = 0;
      }
    } // while
    // close the PDB file
    fin.close();

  } // end of readPQR



  
  vector<string> InPDB::getResSeq() {
    return d_this->resSeq;
  }

  /*void InPDB::changeResSeq(unsigned int i, string newname){
   *  d_this->resSeq[i] = newname;
  }*/


  gmath::Vec InPDB::getAtomPos(unsigned int i) {
    gmath::Vec pos(d_this->X[i], d_this->Y[i], d_this->Z[i]);
    return pos;
  }

  unsigned int InPDB::numAtoms() {
    return d_this->serials.size();
  }

  string InPDB::getResName(unsigned int i) {
    return d_this->resNames[i];
  }

  void InPDB::setResName(int i, string s) {
    d_this->resNames[i] = s;
  }

  unsigned int InPDB::getResNumber(unsigned int i) {
    return d_this->resNums[i];
  }

  string InPDB::getAtomName(unsigned int i) {
    return d_this->atoms[i];
  }
  
  char InPDB::getICode(unsigned int i) {
    return d_this->iCodes[i];
  }
  
  string InPDB::getElement(unsigned int i) {
    return d_this->elements[i];
  }
  
  string InPDB::getCharge(unsigned int i) {
    return d_this->charges[i];
  }

  string InPDB::getType(unsigned int i) {
    return d_this->types[i];
  }

  double InPDB::getOcc(unsigned int i) {
    return d_this->occupancies[i];
  }

  double InPDB::getB(unsigned int i) {
    return d_this->tempFactors[i];
  }

  unsigned int InPDB::getSerial(unsigned int i) {
    return d_this->serials[i];
  }

  char InPDB::getChain(unsigned int i) {
    return d_this->chainIDs[i];
  }
  double InPDB::PQR_getCharges(unsigned int i) {
    return d_this->pqr_charges[i];
  }
  
  double InPDB::PQR_getRadii(unsigned int i) {
    return d_this->pqr_radii[i];
  }

  void InPDB::renumberRes() {
    int oldPDBnum;
    int newPDBnum;
    int resNum = 1;
    for (unsigned int i = 0; i < numAtoms() - 1; ++i) {
      oldPDBnum = d_this->resNums[i];
      newPDBnum = d_this->resNums[i + 1];
      d_this->resNums[i] = resNum;
      if (oldPDBnum != newPDBnum) {
        resNum++;
      }

      //debug:
      //cout << oldPDBnum << "  " << newPDBnum << "  " << resNum << endl;

    }
    d_this->resNums[numAtoms() - 1] = resNum;
  }
  
  
    bool InPDB::hasSerials() {
      return d_this->b_serial;
    }
    
    bool InPDB::hasAtomNames(){
      return d_this->b_atomname;
    }
    
    bool InPDB::hasResNames(){
      return d_this->b_resname;
    }
    
    bool InPDB::hasChainIDs(){
      return d_this->b_chainid;
    }
    
    bool InPDB::hasResSeqs(){
      return d_this->b_resseq;
    }
    
    bool InPDB::hasiCodes(){
      return d_this->b_icode;
    }
    
    bool InPDB::hasX(){
      return d_this->b_x;
    }
    
    bool InPDB::hasY(){
      return d_this->b_y;
    }
    
    bool InPDB::hasZ(){
      return d_this->b_z;
    }
    
    bool InPDB::hasOccupancies(){
      return d_this->b_occupancy;
    }
    
    bool InPDB::hasTempFactors(){
      return d_this->b_tempfactor;
    }
    
    bool InPDB::hasElements(){
      return d_this->b_element;
    }
    
    bool InPDB::hasCharges(){
      return d_this->b_charge;
    }

}

string removeChar(string s, char c) {
  string::iterator it = s.begin();
  for(; it != s.end(); ++it) {
    if (*it == c) {
      s.erase(it);
      --it;
    }
  }
  return s;
}
