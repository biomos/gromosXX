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
#include "Noe.h"

#include <cassert>
#include <cstdlib>
#include <ios>
#include <ostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

#include "VirtualAtom.h"
#include "Neighbours.h"
#include "AtomSpecifier.h"
#include "../gio/StringTokenizer.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/VirtualAtoms.h"
#include "../gromos/Exception.h"

using namespace gio;
using namespace std;
using namespace gcore;
using namespace utils;
using utils::Noe;


class utils::Noe_i{
  friend class utils::Noe;

  const System &d_sys;
  std::vector<VirtualAtom*> d_at[2];
  std::vector<double> d_dist;
  std::vector<double> cor;
  std::vector<int> cortype;
  int d_num;  
  Noe_i(const System &sys): d_sys(sys), d_at(), d_dist(){}
  ~Noe_i(){}
};

//implementation of the new NOE program/constructor, taking a PROADR input
Noe::Noe(System  &sys, const string &line, double dish, double disc):d_this(new Noe_i(sys))
{

  //parse line into Tokens...
  StringTokenizer stok(line, " ");
  std::vector<std::string> tokens = stok.tokenize();
  if (tokens.size() < 13)
    throw Exception("At least 13 input digits are expected! Check input!\n");
  
  // distance
  d_this->d_dist.push_back(atof(tokens[12].c_str()));
  
  // two VirtualAtoms to define an NOE (first with offset 0, second with offset 6)
  int offset = 0;
  for(int k=0; k<2; ++k){
    VirtualAtom::virtual_type type = VirtualAtom::virtual_type(atoi(tokens[4
            +offset].c_str()));

    // int config[4];
    std::vector<int> config(4);
    
    config[0] = atoi(tokens[0+offset].c_str())-1;
    config[1] = atoi(tokens[1+offset].c_str())-1;
    config[2] = atoi(tokens[2+offset].c_str())-1;
    config[3] = atoi(tokens[3+offset].c_str())-1;
    
    //int subtype = atoi(tokens[5+offset].c_str());
    
    d_this->d_at[k].push_back(new VirtualAtom(sys, type, config, dish, disc));
    
    offset = 6;
  }

  // how many virtual atoms per NOE site
  // this is probably deprecated
  d_this->d_num = d_this->d_at[0].size() * d_this->d_at[1].size();
  
}

double Noe::distance(int i)const{
  return distanceVec(i).abs();

}

gmath::Vec Noe::distanceVec(int i)const{
  assert(i<d_this->d_num);
  int at[2];

  at[1] = i % d_this->d_at[1].size();
  at[0] = (i-at[1])/d_this->d_at[1].size();
  
  return d_this->d_at[0][at[0]]->pos() -
	 d_this->d_at[1][at[1]]->pos();

}

const utils::VirtualAtom & Noe::getAtom(int i, int ii) const{
  assert(ii<d_this->d_num);
  int at[2];

  at[1] = ii % d_this->d_at[1].size();
  at[0] = (ii-at[1])/d_this->d_at[1].size();
  
  return *(d_this->d_at[i][at[i]]);
}

string Noe::info(int i)const{
  assert(i<d_this->d_num);
  ostringstream os;

  int at[2];

  at[1] = i % d_this->d_at[1].size();
  at[0] = (i-at[1])/d_this->d_at[1].size();

  for (int j=0 ; j < 2; ++j){
    int atom = d_this->d_at[j][at[j]]->conf().atom(0);
    int mol = d_this->d_at[j][at[j]]->conf().mol(0);
    int type = d_this->d_at[j][at[j]]->type();

    int resNum = d_this->d_sys.mol(mol).topology().resNum(atom);
    string resName = d_this->d_sys.mol(mol).topology().resName(resNum);
    string atName = d_this->d_sys.mol(mol).topology().atom(atom).name();

     
    //std::ostrstream oss;
     ostringstream oss;
    if (type == VirtualAtom::stereo_CH2) {
     oss << " Type 4 Noe! Atoms: "
         << ((d_this->d_at[j][at[j]]->conf().atom(0))+1) << " " 
         << ((d_this->d_at[j][at[j]]->conf().atom(1))+1) << " " 
         << ((d_this->d_at[j][at[j]]->conf().atom(2))+1);
    }

    std::string typefour(oss.str());

    //os.freeze(false);  // Unfreeze so buffer is deleted.

    //cout << typefour << endl;

    switch (type) {
      case VirtualAtom::CH2:
      case VirtualAtom::CH31:
      case VirtualAtom::CH32:
        atName[0] = 'Q';
        break;
      case VirtualAtom::stereo_CH2:
        atName += typefour;
      case VirtualAtom::CH1:
      case VirtualAtom::aromatic:
        atName[0] = 'H';
        break;
    }

    os.setf(ios::right, ios::adjustfield);
    os << setw(5) << resNum+1;
    os.setf(ios::left, ios::adjustfield);
    os <<  setw(5) << resName.c_str() << setw(5) << atName.c_str();
  }

  return string(os.str());
}


int Noe::numDistances()const{
  return d_this->d_num;
}
int Noe::numReferences()const{
  return d_this->d_dist.size();
}


string Noe::distRes(int i)const{

  assert(i<d_this->d_num);
  int at[2];

  at[1] = i % d_this->d_at[1].size();
  at[0] = (i-at[1])/d_this->d_at[1].size();


  ostringstream ss;
  ss.setf(ios::right, ios::adjustfield);
  ss.setf(ios::fixed, ios::floatfield);
  ss.precision(3);

  for (int j = 0; j < 2; ++j) {
    for (int k = 0; k < 4; ++k) {
      int att = d_this->d_at[j][at[j]]->conf().atom(k);
      int mol = d_this->d_at[j][at[j]]->conf().mol(k);
      int offset = 1;
      for (int j = 0; j < mol; ++j)
        offset += d_this->d_sys.mol(j).numAtoms();
      if (att == -1)
        ss << setw(5) << 0;
      else
        ss << setw(5) << att + offset;
    }
    ss << setw(3) << d_this->d_at[j][at[j]]->type() % 7;
  }

  if (i < int(d_this->d_dist.size())) {
    ss << setw(10) << correctedReference(i);
  } else if (i == int(d_this->d_dist.size())) {
    ss << setw(10) << correctedReference(i - 1);
  } else if (i == 2) {
    ss << setw(10) << correctedReference(i - 2);
  } else if (i == 3) {
    ss << setw(10) << correctedReference(i - 3);
  }
  
  ss << setw(10) << 1.0;
  
  return string(ss.str());
}

double Noe::reference(int i)const{
  assert( i < int(d_this->d_dist.size()));
  return d_this->d_dist[i];
}

double Noe::correctedReference(int i)const{
  assert( i < int(d_this->d_dist.size()));
  double cd=d_this->d_dist[i];

  //now come several multiplicity corrections
  //see: Neuhaus D. et.al. J.Bio.NMR, 8 (1996) 292-310

  // for type 5, the experimental distance has to be multiplied by
  // 3^(1/6)
  for(int k=0;k<2;k++){
      if(d_this->d_at[k][0]->type()==VirtualAtom::CH32) {
	cd*=pow(3.0,1.0/6.0); 
      }
  // for type 3, the experimental distance has to be multiplied by
  // 2^(1/6)
      else if(d_this->d_at[k][0]->type()==VirtualAtom::CH2) {
        cd*=pow(2.0,1.0/6.0);
      }
  // for type 6, the experimental distance has to be multiplied by
  // 6^(1/6)
      else if(d_this->d_at[k][0]->type()==VirtualAtom::CH31) {
        cd*=pow(6.0,1.0/6.0);
      }
  }

  //parse the actual read in corrections from the correction file

  for(int k=0;k<2;k++){
    switch(d_this->d_at[k][0]->type()){ 

      case VirtualAtom::CH2:
	cd+=d_this->cor[0];
	break;
      case VirtualAtom::CH32:
	cd+=d_this->cor[1];
	break;
      case VirtualAtom::CH31:
	cd+=d_this->cor[2];
	break;
      default:
	break;
	
    }
  }
  
  return cd;
  
}

double Noe::correction(int type) {

  if ((type < VirtualAtom::CH2) || (type == VirtualAtom::stereo_CH2)){
  ostringstream os;
     os << "GROMOS Noe type not known: "
           << type << "\n";
     throw gromos::Exception("Noe:",os.str());
   }
     
  double t=0;
  for (int i=0; i < int (d_this->cortype.size()); ++i){
    if (d_this->cortype[i] == type) {
      t = d_this->cor[i];}
  }
    
  return t;
}
     

void Noe::setcorrection(int type, double correction) {
  if ((type < VirtualAtom::CH2) || (type == VirtualAtom::stereo_CH2)){
  ostringstream os;
     os << "GROMOS Noe type not known: "
           << type << "\n";
     throw gromos::Exception("Noe:",os.str());
   }
  
  for (int i=0; i < int (d_this->cortype.size()); ++i){
    if (d_this->cortype[i] == type) {
      d_this->cor[i] = correction;}
  }
}

void utils::parse_noelib(const std::vector<std::string>& buffer, std::vector<Noelib>& noelib) {
  if (buffer[0] != "NOELIB")
    throw gromos::Exception("main",
          "NOELIB file does not contain an NOELIB block!");
  if (buffer[buffer.size() - 1].find("END") != 0)
    throw gromos::Exception("prep_noe", "Library file "
          " is corrupted. No END in NOELIB"
          " block. Got\n"
          + buffer[buffer.size() - 1]);

  for (unsigned int j = 1; j < buffer.size() - 1; j++) {
    StringTokenizer tok(buffer[j]);
    vector<string> tokens = tok.tokenize();

    // check for number of tokens and store the data
    if (tokens.size() == 4)
      noelib.push_back(Noelib(tokens[0], tokens[1], tokens[2], tokens[3]));
    else if (tokens.size() == 5) {
      noelib.push_back(Noelib(tokens[0], tokens[1], tokens[2], tokens[3],
              tokens[4]));
    }
  }

  //check for inconsistency in library
  for (int i = 0; i < int (noelib.size()); ++i) {
    const Noelib & A = noelib[i];
    for (int j = 0; j < int (noelib.size()); ++j) {
      const Noelib & B = noelib[j];
      if ((A.resname == B.resname) &&
              (A.orgatomname == B.orgatomname) &&
              (A.gratomname == B.gratomname) &&
              (A.NOETYPE != B.NOETYPE || A.NOESUBTYPE != B.NOESUBTYPE)) {
        std::stringstream sa;
        sa << A.resname << " "
                << A.orgatomname << " "
                << A.gratomname << " "
                << A.NOETYPE;
        if (A.NOESUBTYPE || B.NOESUBTYPE)
          sa << " " << A.NOESUBTYPE;
        sa << " !AND! "
                << B.resname << " "
                << B.orgatomname << " "
                << B.gratomname << " "
                << B.NOETYPE;
        if (A.NOESUBTYPE || B.NOESUBTYPE)
          sa << " " << B.NOESUBTYPE;
        sa << endl;
        throw gromos::Exception("prep_noe ", sa.str() +
                " Inconsistent assigment of NOETYPE within library!");
      }
    }
  }
}

vector<VirtualAtom*> utils::getvirtual(int at, int type, int subtype, System &sys,
        double dish, double disc) {

  int mol = 0, atNum = 0;

  // parse into mol and atom rather than high atom nr.
  while (at >= (atNum += sys.mol(mol).numAtoms())) {
    ++mol;
    if (mol >= sys.numMolecules())
      throw gromos::Exception("prep_noe ", +"Atom number too high in input atom number:\n" + at);
  }
  at -= atNum - sys.mol(mol).numAtoms();

  vector<VirtualAtom*> vat;

  if (type == 4) {
    if (subtype == 0) {

      //we need to automatically generate the r, l atoms...
      vat.push_back(new VirtualAtom(sys, mol, at, VirtualAtom::virtual_type(type), dish, disc));
      vat.push_back(new VirtualAtom(sys, mol, at, VirtualAtom::virtual_type(type), dish, disc, 1));
    } else if (subtype == 1) {
      vat.push_back(new VirtualAtom(sys, mol, at, VirtualAtom::virtual_type(type), dish, disc));
    } else if (subtype == 2) {
      vat.push_back(new VirtualAtom(sys, mol, at, VirtualAtom::virtual_type(type), dish, disc, 1));
    }
  } else if (type == -1 || type == -2) {
    vat.push_back(new VirtualAtom("noe", sys, mol, at, VirtualAtom::virtual_type(type), subtype, dish, disc));
  } else vat.push_back(new VirtualAtom(sys, mol, at, VirtualAtom::virtual_type(type), dish, disc));

  return vat;

}
