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
#include "AminoAcid.h"

#include <string>
#include <iostream>
#include <iomanip>
#include <map>
#include <utility>
#include <vector>
#include <sstream>

#include "../gromos/Exception.h"

using namespace std;

namespace utils {

  void gromosAminoAcidLibrary::loadHardcoded45A4(void) {

    // clear the library if not empty
    if(lib.size() > 0) {
      lib.clear();
    }

    string pdbname;
    vector<string> donors;
    vector<string> acceptors;
    gromosAminoAcid gaa;

    // set the library version
    version = "45A4";

    pdbname = "ALA";
    gaa.acid = "ALA";
    gaa.base = "ALA";
    gaa.pKa = 2.33;
    gaa.pKb = 9.71;
    gaa.pKc = -1.0;
    donors.push_back("N");
    gaa.Hdonors.insert(pair<string, vector<string> > ("ALA", donors));
    acceptors.push_back("O");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("ALA", acceptors));
    donors.clear();
    acceptors.clear();
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "ARG";
    gaa.acid = "ARG";
    gaa.base = "ARGN";
    gaa.pKa = 2.03;
    gaa.pKb = 9.00;
    gaa.pKc = 12.10;
    donors.push_back("N");
    donors.push_back("NE");
    donors.push_back("NH1");
    donors.push_back("NH2");
    gaa.Hdonors.insert(pair<string, vector<string> > ("ARG", donors));
    acceptors.push_back("O");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("ARG", acceptors));
    donors.clear();
    acceptors.clear();
    donors.push_back("N");
    donors.push_back("NE");
    donors.push_back("NH1");
    donors.push_back("NH2");
    gaa.Hdonors.insert(pair<string, vector<string> > ("ARGN", donors));
    acceptors.push_back("O");
    acceptors.push_back("NH1"); // check again later
    gaa.Hacceptors.insert(pair<string, vector<string> > ("ARGN", acceptors));
    donors.clear();
    acceptors.clear();
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "ASN";
    gaa.acid = "ASN";
    gaa.base = "ASN";
    gaa.pKa = 2.16;
    gaa.pKb = 8.73;
    gaa.pKc = -1.0;
    donors.push_back("N");
    donors.push_back("ND2");
    gaa.Hdonors.insert(pair<string, vector<string> > ("ASN", donors));
    acceptors.push_back("O");
    acceptors.push_back("OD1");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("ASN", acceptors));
    donors.clear();
    acceptors.clear();
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "ASP";
    gaa.acid = "ASPH";
    gaa.base = "ASP";
    gaa.pKa = 1.95;
    gaa.pKb = 9.66;
    gaa.pKc = 3.71;
    donors.push_back("N");
    gaa.Hdonors.insert(pair<string, vector<string> > ("ASP", donors));
    acceptors.push_back("O");
    acceptors.push_back("OD1");
    acceptors.push_back("OD2");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("ASP", acceptors));
    donors.clear();
    acceptors.clear();
    donors.push_back("N");
    donors.push_back("OD2");
    gaa.Hdonors.insert(pair<string, vector<string> > ("ASPH", donors));
    acceptors.push_back("O");
    acceptors.push_back("OD2"); // check later again
    acceptors.push_back("OD1");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("ASPH", acceptors));
    donors.clear();
    acceptors.clear();
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "CYS";
    gaa.acid = "CYSH";
    gaa.base = "CYS";
    gaa.pKa = 1.91;
    gaa.pKb = 10.28;
    gaa.pKc = 8.14;
    donors.push_back("N");
    gaa.Hdonors.insert(pair<string, vector<string> > ("CYS", donors));
    acceptors.push_back("O");
    acceptors.push_back("SG");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("CYS", acceptors));
    donors.clear();
    acceptors.clear();
    donors.push_back("N");
    donors.push_back("SG");
    gaa.Hdonors.insert(pair<string, vector<string> > ("CYSH", donors));
    acceptors.push_back("O");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("CYSH", acceptors));
    donors.clear();
    acceptors.clear();
    donors.push_back("N");
    gaa.Hdonors.insert(pair<string, vector<string> > ("CYS1", donors));
    acceptors.push_back("O");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("CYS1", acceptors));
    donors.clear();
    acceptors.clear();
    donors.push_back("N");
    gaa.Hdonors.insert(pair<string, vector<string> > ("CYS2", donors));
    acceptors.push_back("O");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("CYS2", acceptors));
    donors.clear();
    acceptors.clear();
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "GLN";
    gaa.acid = "GLN";
    gaa.base = "GLN";
    gaa.pKa = 2.18;
    gaa.pKb = 9.00;
    gaa.pKc = -1.0;
    donors.push_back("N");
    donors.push_back("NE2");
    gaa.Hdonors.insert(pair<string, vector<string> > ("GLN", donors));
    acceptors.push_back("O");
    acceptors.push_back("OE");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("GLN", acceptors));
    donors.clear();
    acceptors.clear();
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "GLU";
    gaa.acid = "GLUH";
    gaa.base = "GLU";
    gaa.pKa = 2.16;
    gaa.pKb = 9.58;
    gaa.pKc = 4.15;
    donors.push_back("N");
    gaa.Hdonors.insert(pair<string, vector<string> > ("GLU", donors));
    acceptors.push_back("O");
    acceptors.push_back("OE1");
    acceptors.push_back("OE2");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("GLU", acceptors));
    donors.clear();
    acceptors.clear();
    donors.push_back("N");
    donors.push_back("OE2");
    gaa.Hdonors.insert(pair<string, vector<string> > ("GLUH", donors));
    acceptors.push_back("O");
    acceptors.push_back("OE2"); // check later again
    acceptors.push_back("OE1");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("GLUH", acceptors));
    donors.clear();
    acceptors.clear();
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "GLY";
    gaa.acid = "GLY";
    gaa.base = "GLY";
    gaa.pKa = 2.34;
    gaa.pKb = 9.58;
    gaa.pKc = -1.0;
    donors.push_back("N");
    gaa.Hdonors.insert(pair<string, vector<string> > ("GLY", donors));
    acceptors.push_back("O");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("GLY", acceptors));
    donors.clear();
    acceptors.clear();
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "HIS";
    gaa.acid = "HISH";
    gaa.base = "HISX";
    gaa.pKa = 1.70;
    gaa.pKb = 9.09;
    gaa.pKc = 6.04;
    donors.push_back("N");
    donors.push_back("ND1");
    donors.push_back("NE2");
    gaa.Hdonors.insert(pair<string, vector<string> > ("HISH", donors));
    acceptors.push_back("O");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("HISH", acceptors));
    donors.clear();
    acceptors.clear();
    donors.push_back("N");
    donors.push_back("ND1");
    gaa.Hdonors.insert(pair<string, vector<string> > ("HISA", donors));
    acceptors.push_back("O");
    acceptors.push_back("NE2");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("HISA", acceptors));
    donors.clear();
    acceptors.clear();
    donors.push_back("N");
    donors.push_back("NE2");
    gaa.Hdonors.insert(pair<string, vector<string> > ("HISB", donors));
    acceptors.push_back("O");
    acceptors.push_back("ND1");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("HISB", acceptors));
    donors.clear();
    acceptors.clear();
    donors.push_back("N");
    gaa.Hdonors.insert(pair<string, vector<string> > ("HISX", donors)); // come back here and have fun ;-)
    acceptors.push_back("O");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("HISX", acceptors));
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    donors.clear();
    acceptors.clear();
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "ILE";
    gaa.acid = "ILE";
    gaa.base = "ILE";
    gaa.pKa = 2.26;
    gaa.pKb = 9.60;
    gaa.pKc = -1.0;
    donors.push_back("N");
    gaa.Hdonors.insert(pair<string, vector<string> > ("ILE", donors));
    acceptors.push_back("O");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("ILE", acceptors));
    donors.clear();
    acceptors.clear();
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "LEU";
    gaa.acid = "LEU";
    gaa.base = "LEU";
    gaa.pKa = 2.32;
    gaa.pKb = 9.58;
    gaa.pKc = -1.0;
    donors.push_back("N");
    gaa.Hdonors.insert(pair<string, vector<string> > ("LEU", donors));
    acceptors.push_back("O");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("LEU", acceptors));
    donors.clear();
    acceptors.clear();
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "LYS";
    gaa.acid = "LYSH";
    gaa.base = "LYS";
    gaa.pKa = 2.15;
    gaa.pKb = 9.16;
    gaa.pKc = 10.67;
    donors.push_back("N");
    donors.push_back("NZ");
    gaa.Hdonors.insert(pair<string, vector<string> > ("LYS", donors));
    acceptors.push_back("O");
    acceptors.push_back("NZ"); // check later again
    gaa.Hacceptors.insert(pair<string, vector<string> > ("LYS", acceptors));
    donors.clear();
    acceptors.clear();
    donors.push_back("N");
    donors.push_back("NZ");
    gaa.Hdonors.insert(pair<string, vector<string> > ("LYSH", donors));
    acceptors.push_back("O");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("LYSH", acceptors));
    donors.clear();
    acceptors.clear();
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "MET";
    gaa.acid = "MET";
    gaa.base = "MET";
    gaa.pKa = 2.16;
    gaa.pKb = 9.08;
    gaa.pKc = -1.0;
    donors.push_back("N");
    gaa.Hdonors.insert(pair<string, vector<string> > ("MET", donors));
    acceptors.push_back("O");
    acceptors.push_back("SD");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("MET", acceptors));
    donors.clear();
    acceptors.clear();
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "PHE";
    gaa.acid = "PHE";
    gaa.base = "PHE";
    gaa.pKa = 2.18;
    gaa.pKb = 9.09;
    gaa.pKc = -1.0;
    donors.push_back("N");
    gaa.Hdonors.insert(pair<string, vector<string> > ("PHE", donors));
    acceptors.push_back("O");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("PHE", acceptors));
    donors.clear();
    acceptors.clear();
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "PRO";
    gaa.acid = "PRO";
    gaa.base = "PRO";
    gaa.pKa = 1.95;
    gaa.pKb = 10.47;
    gaa.pKc = -1.0;
    // no donors for PRO! => add the empty donor vector donors
    gaa.Hdonors.insert(pair<string, vector<string> > ("PRO", donors));
    acceptors.push_back("O");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("PRO", acceptors));
    acceptors.clear();
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "SER";
    gaa.acid = "SER";
    gaa.base = "SER";
    gaa.pKa = 2.13;
    gaa.pKb = 9.05;
    gaa.pKc = -1.0;
    donors.push_back("N");
    donors.push_back("OG");
    gaa.Hdonors.insert(pair<string, vector<string> > ("SER", donors));
    acceptors.push_back("O");
    acceptors.push_back("OH"); // think about later
    gaa.Hacceptors.insert(pair<string, vector<string> > ("SER", acceptors));
    donors.clear();
    acceptors.clear();
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "THR";
    gaa.acid = "THR";
    gaa.base = "THR";
    gaa.pKa = 2.20;
    gaa.pKb = 8.96;
    gaa.pKc = -1.0;
    donors.push_back("N");
    donors.push_back("OG1");
    gaa.Hdonors.insert(pair<string, vector<string> > ("THR", donors));
    acceptors.push_back("O");
    acceptors.push_back("OG1"); // think about later
    gaa.Hacceptors.insert(pair<string, vector<string> > ("THR", acceptors));
    donors.clear();
    acceptors.clear();
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "TRP";
    gaa.acid = "TRP";
    gaa.base = "TRP";
    gaa.pKa = 2.38;
    gaa.pKb = 9.34;
    gaa.pKc = -1.0;
    donors.push_back("N");
    donors.push_back("NH1");
    gaa.Hdonors.insert(pair<string, vector<string> > ("TRP", donors));
    acceptors.push_back("O");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("TRP", acceptors));
    donors.clear();
    acceptors.clear();
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "TYR";
    gaa.acid = "TYR";
    gaa.base = "TYR";
    gaa.pKa = 2.24;
    gaa.pKb = 9.04;
    gaa.pKc = -1.0;
    donors.push_back("N");
    donors.push_back("OH");
    gaa.Hdonors.insert(pair<string, vector<string> > ("TYR", donors));
    acceptors.push_back("O");
    acceptors.push_back("OH"); // think about later
    gaa.Hacceptors.insert(pair<string, vector<string> > ("TYR", acceptors));
    donors.clear();
    acceptors.clear();
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "VAL";
    gaa.acid = "VAL";
    gaa.base = "VAL";
    gaa.pKa = 2.27;
    gaa.pKb = 9.52;
    gaa.pKc = -1.0;
    donors.push_back("N");
    gaa.Hdonors.insert(pair<string, vector<string> > ("VAL", donors));
    acceptors.push_back("O");
    gaa.Hacceptors.insert(pair<string, vector<string> > ("VAL", acceptors));
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    
  }


  void gromosAminoAcidLibrary::writeLibrary(ostream &os, string title) {

    // only wirte a library if there is something to write
    if (lib.size() > 0) {

      os << "TITLE\n";
      os << title << endl;
      os << "END\n";
      os << "VERSION\n";
      os << version << endl;
      os << "END\n";
      for (map<string, gromosAminoAcid>::iterator it = lib.begin();
              it != lib.end(); it++) {
        os << "AMINOACID\n";
        os << "#" << setw(11) << "PDB name" << setw(12) << "acid name" << setw(12)
                << "base name" << endl;
        os << setw(12) << it->first
                << setw(12) << it->second.acid
                << setw(12) << it->second.base << endl;
        os << "#\n";
        os << "#" << setw(11) << "pKa" << setw(12) << "pKb" << setw(12) << "pKc" << endl;
        os << setw(12) << it->second.pKa << setw(12) << it->second.pKb
                << setw(12) << it->second.pKc << endl;
        /*
        os << "# number of H-donors" << endl;
        os << setw(12) << it->second.Hdonors.size() << endl;
        os << "#" << setw(11) << "residue" << setw(12) << "H-donors" << endl;
        multimap<string, string>::iterator iit = it->second.Hdonors.begin();
        string resName = iit->first;
        os << setw(12) << resName;
        while(iit != it->second.Hdonors.end()) {
          if(resName == iit->first) {
            os << setw(12) << iit->second;
          } else {
            resName = iit->first;
            os << endl << setw(12) << resName << setw(12) << iit->second;
          }
          iit++;
        }
        os << endl;
        os << "# number of H-acceptors" << endl;
        os << setw(12) << it->second.Hacceptors.size() << endl;
        os << "#" << setw(11) << "residue" << setw(12) << "H-acceptors" << endl;
        */
        os << "END\n#\n";
      }

    } else {

    }

  }

  string gromosAminoAcidLibrary::pdb2acid(std::string PDBname) {
    string s = "XXX";
    if(lib.find(PDBname) != lib.end()) {
      s = lib.find(PDBname)->second.acid;
    }
    return s;
  }

  string gromosAminoAcidLibrary::pdb2base(std::string PDBname) {
    string s = "XXX";
    if(lib.find(PDBname) != lib.end()) {
      s = lib.find(PDBname)->second.base;
    }
    return s;
  }

  double gromosAminoAcidLibrary::pKa(string PDBname) {
    double value = -1;
    if(lib.find(PDBname) != lib.end()) {
      value = lib.find(PDBname)->second.pKa;
    }
    return value;
  }

  double gromosAminoAcidLibrary::pKb(string PDBname) {
    double value = -1;
    if(lib.find(PDBname) != lib.end()) {
      value = lib.find(PDBname)->second.pKb;
    }
    return value;
  }

  double gromosAminoAcidLibrary::pKc(string PDBname) {
    double value = -1;
    if(lib.find(PDBname) != lib.end()) {
      value = lib.find(PDBname)->second.pKc;
    }
    return value;
  }

  vector<string> gromosAminoAcidLibrary::rHdonors(std::string PDBname, std::string GROMOSname) {
    vector<string> s;

    if(lib.find(PDBname) != lib.end()) {
      if (lib.find(PDBname)->second.Hdonors.find(GROMOSname) != lib.find(PDBname)->second.Hdonors.end()) {
      s = lib.find(PDBname)->second.Hdonors.find(GROMOSname)->second;
      } else{
        stringstream msg;
        msg << "no entry for the acceptors of the amino acid " << GROMOSname << " (GROMOS name) in the amino acid library";
        throw gromos::Exception("AminoAcid", msg.str());
      }
    }else{
      stringstream msg;
      msg << "no entry for the amino acid " << PDBname << " (PDB name) in the amino acid library";
      throw gromos::Exception("AminoAcid", msg.str());
    }
    return s;
  }

  vector<string> gromosAminoAcidLibrary::rHacceptors(std::string PDBname, std::string GROMOSname) {
    vector<string> s;

    if (lib.find(PDBname) != lib.end()) {
      if (lib.find(PDBname)->second.Hacceptors.find(GROMOSname) != lib.find(PDBname)->second.Hacceptors.end()) {
        s = lib.find(PDBname)->second.Hacceptors.find(GROMOSname)->second;
      } else {
        stringstream msg;
        msg << "no entry for the donors of the amino acid " << GROMOSname << " (GROMOS name) in the amino acid library";
        throw gromos::Exception("AminoAcid", msg.str());
      }
    } else {
      stringstream msg;
      msg << "no entry for the amino acid " << PDBname << " (PDB name) in the amino acid library";
      throw gromos::Exception("AminoAcid", msg.str());
    }
    
    return s;
  }
  
  std::map<std::string, gromosAminoAcid> gromosAminoAcidLibrary::getAminoAcids() {
    return lib;
  }

}
