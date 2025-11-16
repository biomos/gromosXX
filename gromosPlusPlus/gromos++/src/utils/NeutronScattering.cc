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
#include "NeutronScattering.h"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <ios>
#include <ostream>
#include <string>
#include <utility>
#include <vector>
#include <map>
#include <set>
#include <iomanip>
#include <sstream>
#include <fstream>

#include "RDF.h"
#include "AtomSpecifier.h"
#include "../gromos/Exception.h"
#include "../args/Arguments.h"
#include "../gcore/System.h"
#include "../gmath/Physics.h"
#include "../gio/Ginstream.h"

using namespace std;
using namespace args;
using namespace gcore;
using namespace gmath;

namespace utils {

  class iNS {
  public:

    int d_grid;
    double d_cut;
    double d_Qmax;
    vector<RDF> d_rdf;
    vector<vector<double> > d_Sintra;
    vector<vector<double> > d_Sinter;
    vector<double> d_intensity;
    vector<double> IofQ;
    System *d_sys;
    multimap<int, int> d_comb;
    vector<double> d_intraC2Wdist;
    map<int, double> d_scattLen;
    map<int, map<int, double> > d_sigma;
    AtomSpecifier d_atoms;
    vector<double> d_weightIntra;
    vector<double> d_weightInter;
    map<int, double> d_afraction;
    const args::Arguments *d_args;

  };

  NS::NS(System *sys, const args::Arguments *args) {
    d_this = new iNS;
    setSystem(sys);
    d_this->d_grid = 200;
    d_this->d_cut = 1.5;
    d_this->d_Qmax = 200;
    d_this->d_args = args;
  }

  NS::~NS(void) {
    if (d_this) {
      delete d_this;
    }
  }

  void NS::setGrid(int grid) {
    assert(d_this != NULL);
    assert(d_this->d_rdf.size() == d_this->d_Sinter.size() &&
            d_this->d_rdf.size() == d_this->d_Sintra.size());
    d_this->d_grid = grid;
    for (unsigned int i = 0; i < d_this->d_rdf.size(); ++i) {
      d_this->d_rdf[i].setGrid(grid);
      d_this->d_Sinter[i].resize(grid);
      d_this->d_Sintra[i].resize(grid);
    }
    d_this->d_intensity.resize(grid);
  }

  void NS::setCut(double cut) {
    assert(d_this != NULL);
    d_this->d_cut = cut;
    for (unsigned int i = 0; i < d_this->d_rdf.size(); ++i) {
      d_this->d_rdf[i].setCut(cut);
    }
  }

  void NS::setQmax(int Qmax) {
    assert(d_this != NULL);
    d_this->d_Qmax = Qmax;
  }

  int NS::addAtoms(string s) {
    assert(d_this != NULL);
    d_this->d_atoms.addSpecifier(s);
    d_this->d_atoms.sort();
    return d_this->d_atoms.size();
  }

  void NS::setSystem(System *sys) {
    assert(d_this != NULL);
    d_this->d_sys = sys;
    d_this->d_atoms.setSystem(*sys);
    for(unsigned int i = 0; i < d_this->d_rdf.size(); ++i) {
      d_this->d_rdf[i].setSystem(sys);
    }
  }

  int NS::getCombinations(void) {
    // make sure there are no (old) combinations in the list now
    d_this->d_comb.clear();
    for(unsigned int c = 0; c < d_this->d_atoms.size(); ++c) {
      for(unsigned int w = 0; w < d_this->d_atoms.size(); ++w) {
        int iacc = d_this->d_atoms.iac(c); // has to be in this loop since it
                                           // need to be reset in case of
                                           // iacw > iacc, which is switched
                                           // below
        int iacw = d_this->d_atoms.iac(w);
        // to make sure the centres are smaller than the with in the d_comb map...
        if(iacw < iacc) {
          int tmp = iacc;
          iacc = iacw;
          iacw = tmp;
        }
        bool found = false;
        multimap<int, int>::const_iterator start = d_this->d_comb.lower_bound(iacc);
        multimap<int, int>::const_iterator stop = d_this->d_comb.upper_bound(iacc);
        multimap<int, int>::const_iterator it;
        for (it = start; it != stop; ++it) {
          if (it->second == iacw) {
            found = true;
            break;
          }
        }
        if (!found) {
          d_this->d_comb.insert(pair<int, int>(iacc, iacw));
        }
      }
    }
    // sort the multimap of combinations
    map<int, set<int> > comb;
    {
      multimap<int, int>::iterator it;
      for(it = d_this->d_comb.begin(); it != d_this->d_comb.end(); ++it) {
        map<int, set<int> >::iterator iit = comb.find(it->first);
        if(iit != comb.end()) {
          iit->second.insert(it->second);
        } else {
          set<int> s;
          s.insert(it->second);
          comb.insert(pair<int, set<int> >(it->first, s));
        }
      }
      d_this->d_comb.clear();
      map<int, set<int> >::iterator itc;
      for(itc = comb.begin(); itc != comb.end(); itc++) {
        set<int>::iterator its;
        for(its = itc->second.begin(); its != itc->second.end(); ++its) {
          d_this->d_comb.insert(pair<int, int>(itc->first, *its));
        }
      }
    }
    // now resize the depending vector lengths
    d_this->d_Sinter.resize(d_this->d_comb.size());
    d_this->d_Sintra.resize(d_this->d_comb.size());
    d_this->d_rdf.resize(d_this->d_comb.size());
    d_this->d_weightInter.resize(d_this->d_comb.size());
    d_this->d_weightIntra.resize(d_this->d_comb.size());
    d_this->d_intraC2Wdist.resize(d_this->d_comb.size());
    // now calculate the mole fractions
    d_this->d_afraction.clear();
    map<int, set<int> >::iterator it;
    for(it = comb.begin(); it != comb.end(); ++it) {
      int num = 0;
      for(unsigned int a = 0; a < d_this->d_atoms.size(); ++a) {
        if(it->first == d_this->d_atoms.iac(a)) {
          num++;
        }
      }
      d_this->d_afraction.insert(pair<int, double>(it->first, (double)num/(double)d_this->d_atoms.size()));
    }
    return d_this->d_comb.size();
  }

  void NS::check(void) {
    if(!d_this) {
      stringstream msg;
      msg << "inertialisation of the implementation calss iNS failed";
      throw gromos::Exception("class utils::NS", msg.str());
    }
    if(!d_this->d_sys) {
      stringstream msg;
      msg << "no system (gcore::System) set";
      throw gromos::Exception("class utils::NS", msg.str());
    }
  }

  void NS::setRDFatoms() {
    assert(d_this != NULL);
    assert(d_this->d_comb.size() == d_this->d_rdf.size() && d_this->d_comb.size() > 0);
    int i = 0;
    multimap<int, int>::iterator it;
    for(it = d_this->d_comb.begin(); it != d_this->d_comb.end(); ++it) {
      for(unsigned int asn = 0; asn < d_this->d_atoms.size(); ++asn) {
        if (d_this->d_atoms.iac(asn) == it->first) {
          int m = d_this->d_atoms.mol(asn);
          int a = d_this->d_atoms.atom(asn);
          d_this->d_rdf[i].addCentersAtom(m, a);
        }
        if (d_this->d_atoms.iac(asn) == it->second) {
          int m = d_this->d_atoms.mol(asn);
          int a = d_this->d_atoms.atom(asn);
          d_this->d_rdf[i].addWithAtom(m, a);
        }
      }
      ++i;
    }
  }

  void NS::calcRDFsInterAll() {
    assert(d_this != NULL);
    assert(d_this->d_comb.size() == d_this->d_rdf.size());
    for(unsigned int i = 0; i < d_this->d_rdf.size(); ++i) {
      d_this->d_intraC2Wdist[i] = d_this->d_rdf[i].calculateInterPDens(d_this->d_atoms.size());
    }
  }

  void NS::printRDFs(std::ostream &os) {
    for(unsigned int i = 0; i < d_this->d_rdf.size(); ++i) {
      d_this->d_rdf[i].print(os);
      os << endl;
    }
  }

  void NS::printIntensity(std::ostream& os) {
    double Qmin = 2 * physConst.get_pi() * double(d_this->d_grid) /
            ((double(d_this->d_grid) - 0.5) * d_this->d_cut);
    double dQ = (d_this->d_Qmax - Qmin) / d_this->d_grid;
    os.precision(9);
    for(int g = 0; g < d_this->d_grid; ++g) {
      os << scientific << setw(20) << Qmin + g * dQ << setw(20)
              << d_this->d_intensity[g] << endl;
    }
  }

  void NS::getWeights(void) {
    assert(d_this != NULL);
    assert(d_this->d_atoms.size() != 0);
    unsigned int i = 0;
    multimap<int, int>::iterator it;
    for (it = d_this->d_comb.begin(); it != d_this->d_comb.end(); ++it) {
      // make sure there are no old weights
      d_this->d_weightIntra[i] = 0;
      d_this->d_weightInter[i] = 0;
      // get the intramolecular weights
      for (unsigned int c = 0; c < d_this->d_atoms.size(); ++c) {
        for (unsigned int w = 0; w < d_this->d_atoms.size(); ++w) {
          int iacc = d_this->d_atoms.iac(c); // has to be in this loop since it
                                             // need to be reset in case of
                                             // iacw > iacc, which is switched
                                             // below
          int iacw = d_this->d_atoms.iac(w);
          // make sure iacc <= iacw
          if (iacw < iacc) {
            int tmp = iacc;
            iacc = iacw;
            iacw = tmp;
          }
          // skip if centre and with atom are in different molecules
          if (d_this->d_atoms.mol(c) != d_this->d_atoms.mol(w)) {
            continue;
          }
          // skip if the same molecule but also the same atom
          if(d_this->d_atoms.mol(c) == d_this->d_atoms.mol(w) &&
                  d_this->d_atoms.atom(c) == d_this->d_atoms.atom(w)) {
            continue;
          }
          // skipt if centre or with atoms have not the right type (IAC)
          if(iacc != it->first || iacw != it->second) {
            continue;
          }
          d_this->d_weightIntra[i]++;
        }
      }
      // divide by the total number of atoms (we want the intensity per atom)
      d_this->d_weightIntra[i] /= d_this->d_atoms.size();
      // and now the inter-molecular weights (we only loop over the atom types
      int f = 0;
      if (it->first == it->second) {
        f = 1;
      }
      double afraci = d_this->d_afraction.find(it->first)->second;
      double afracj = d_this->d_afraction.find(it->second)->second;
      d_this->d_weightInter[i] = (2.0 - f) * afraci * afracj;
      i++;
    }
  }

  /**
   * Prints the combination of centre to with IAC numbers.
   */
  void NS::print(ostream &os) {

    os << "NEUTRONSCATTERING\n";
    os << "# IACI: IAC number of atom i\n";
    os << "# IACJ: IAC number of atom j\n";
    os << "# WIJM: weight of the corresponding intra-molecular structure factor\n";
    os << "# WIJD: weight of the corresponding inter-molecular structure factor\n";
    os << "# NC  : number of IAC combinations (centre-to-with atoms)\n";
    os << "#\n";
    os << "#" << setw(7) << "NC" << endl;
    os << setw(8) << d_this->d_comb.size() << endl;
    os << "#" << setw(7) << "IACI" << setw(8) << "IACJ"
            << setw(20) << "WIJM" << setw(20) << "WIJD" << endl;
    os.precision(9);
    int i = 0;
    for (multimap<int, int>::iterator it = d_this->d_comb.begin();
            it != d_this->d_comb.end(); ++it) {
      os << setw(8) << (it->first) + 1 << setw(8) << (it->second) + 1<< scientific
              << setw(20) << d_this->d_weightIntra[i]
              << setw(20) << d_this->d_weightInter[i] << endl;
      ++i;
    }
    os << "# NDAT: Number of different atom types\n";
    os << "# IAC:  Integer atom code (atom tupe)\n";
    os << "# SCL:  Scattering lenghts for atoms of type IAC\n";
    os << "#\n";
    os << "#" << setw(14) << "NDAT" << setw(15) << endl;
    os << setw(15) << d_this->d_afraction.size() << endl;
    os << "#" << setw(14) << "IAC" << setw(15) << setw(20) << "SCL"
            << setw(20) << "AFR" << endl;
    for(map<int, double>::iterator it = d_this->d_afraction.begin();
            it != d_this->d_afraction.end(); ++it) {
      double sl;
      if(d_this->d_scattLen.find(it->first) != d_this->d_scattLen.end()) {
        sl = d_this->d_scattLen.find(it->first)->second;
      } else {
        stringstream msg;
        msg << "no scattering length defined for IAC = " << it->first + 1;
        throw gromos::Exception("NeutronScattering", msg.str());
      }
      os << setw(15) << (it->first) + 1
              << setw(20) << scientific
              << sl
              << setw(20) << it->second << endl;
    }
    os << "END\n";
  }

  void NS::readScattlen(string fname) {
    string s;
    stringstream ss;
    bool inSL = false; // to check if we are in the block or not
    // open the file for reading
    ifstream fin(fname.c_str());
    // read line by line
    while (true) {
      getline(fin, s);
      // error message if no block SCATTLENGTHS or no END in block found
      if (fin.eof()) {
        if (inSL) {
          stringstream msg;
          msg << "no END in SCATTLENGTHS block of file " << fname;
          throw gromos::Exception("NeutronScattering.cc", msg.str());
        } else {
          stringstream msg;
          msg << "no SCATTLENGTHS block in file " << fname;
          throw gromos::Exception("NeutronScattering.cc", msg.str());
        }
      }
      // get rid of comments
      s = s.substr(0, s.find('#'));
      // finish if we reach the END of the SCATTLENGTHS block
      if (s == "END" && inSL) {
        break;
      }
      // add the line if there is something left and we are in SCATTLENGTHS
      if (s != "" && inSL) {
        ss << s << endl;
      }
      // from now on we are in SCATTLENGTHS
      if (s == "SCATTLENGTHS") {
        inSL = true;
      }
    }
    fin.close();
    // make sure there are no old scattering lengths saved
    d_this->d_scattLen.clear();
    // now save the read data
    int iac;
    double sl;
    stringstream inp;
    while (ss >> s) {
      inp << s;
      inp >> iac;
      iac--;
      if (inp.fail() || inp.bad()) {
        stringstream msg;
        msg << "expected a positive integer number when reading the IAC but found " << s;
        throw gromos::Exception("NeutronScattering.cc", msg.str());
      }
      if(!(ss >> s)) {
        stringstream msg;
        msg << "bad line in SCATTLENGTHS block, could not read scattering length for IAC " << iac + 1;
        throw gromos::Exception("NeutronScattering.cc", msg.str());
      }
      inp.clear();
      inp.str("");
      inp << s;
      inp >> sl;
      if (inp.fail() || inp.bad()) {
        stringstream msg;
        msg << "scattering lengths " << s << " for IAC " << iac + 1 << " cannot be converted into a double number";
        throw gromos::Exception("NeutronScattering.cc", msg.str());
      }
      inp.clear();
      inp.str("");
      // now add the scattering length to the map
      if(d_this->d_scattLen.find(iac) == d_this->d_scattLen.end()) {
        d_this->d_scattLen.insert(pair<int, double>(iac, sl));
      } else if(d_this->d_scattLen.find(iac)->second != sl) {
        stringstream msg;
        msg << "two different scattering lengths defined in " << fname << " for atom type with IAC = " << iac + 1;
        throw gromos::Exception("NeutronScattering.cc", msg.str());
      }
    }
  }

  void NS::readSigma(string fname) {

    // make sure there is no old stuff
    d_this->d_sigma.clear();

    // read file into a buffer
    gio::Ginstream file(fname);
    vector<string> buffer;
    file.getblock(buffer);
    file.close();

    // check if the SIGMA block is there and complett
    if (buffer[0] != "SIGMAS") {
      throw gromos::Exception("NeutronScattering", "Could not read SIGMAS block in" +
              file.name());
    }
    if (buffer[buffer.size() - 1].find("END") != 0) {
      throw gromos::Exception("NeutronScattering", file.name() + "is corrupt: no END"
              " in " + buffer[0] + " block. Got \n" + buffer[buffer.size() - 1]);
    }

    vector<string>::iterator iter = buffer.begin() + 1;
    istringstream is(*iter);

    if (buffer.size() <= 2) {
      throw gromos::Exception("NeutronScattering", file.name() +
              " gives no information");
    }

    vector<double> tmp(3); // iaci, iacj, sigma

    for (unsigned int i = 0; i < buffer.size() - 2; i++, ++iter) {
      is.clear();
      is.str(*iter);
      if (!(is >> tmp[0] >> tmp[1] >> tmp[2])) {
        throw gromos::Exception("NeutronScattering", "Error reading file" + file.name());
      }
      tmp[0]--;
      tmp[1]--;
      map<int, double> t;
      t.insert(pair<int, double>(int(tmp[1]), tmp[2]));
      if(d_this->d_sigma.find(int(tmp[0])) != d_this->d_sigma.end()) {
        d_this->d_sigma.find(int(tmp[0]))->second.insert(pair<int, double>(int(tmp[1]), tmp[2]));
      } else {
        d_this->d_sigma.insert(pair<int, map<int, double> >(int(tmp[0]), t));
      }
    }
  }

  /* end of NeutronScattering::read_sigma(string fname) */

  double NS::sinc(double x) {
    return sin(x) / x;
  }

  void NS::calcSintra(void) {
    double Qmin = 2 * physConst.get_pi() * double(d_this->d_grid) /
            ((double(d_this->d_grid) - 0.5) * d_this->d_cut);
    double dQ = (d_this->d_Qmax - Qmin) / double(d_this->d_grid - 1);
    multimap<int, int>::iterator it;
    int count = 0;
    for (it = d_this->d_comb.begin(); it != d_this->d_comb.end(); ++it) {
      double r_ij = d_this->d_intraC2Wdist[count];
      if(d_this->d_sigma.find(it->first) == d_this->d_sigma.end() ||
              d_this->d_sigma.find(it->first)->second.find(it->second) == d_this->d_sigma.find(it->first)->second.end()) {
        stringstream msg;
        msg << "no positional root-mean-square deviation (sigma) found for atom IAC combination " << it->first << " " << it->second << endl;
        throw gromos::Exception("NeutronScattering.cc", msg.str());
      }
      double sigma = d_this->d_sigma.find(it->first)->second.find(it->second)->second;
      if (r_ij > 0) {
        for (int g = 0; g < d_this->d_grid; ++g) {
          double arg = sigma * (Qmin + g * dQ);
          d_this->d_Sintra[count][g] = sinc(d_this->d_intraC2Wdist[count] *
                  (Qmin + g * dQ)) * exp(-(arg * arg) / 2);
        }
      }
      count++;
    }
  }

  void NS::calcSintraElastic(void) {
    double Qmin = 2 * physConst.get_pi() * double(d_this->d_grid) /
            ((double(d_this->d_grid) - 0.5) * d_this->d_cut);
    double dQ = (d_this->d_Qmax - Qmin) / double(d_this->d_grid - 1);
    for (unsigned int i = 0; i < d_this->d_rdf.size(); ++i) {
      double r_ij = d_this->d_intraC2Wdist[i];
      if (r_ij > 0) {
        for (int g = 0; g < d_this->d_grid; ++g) {
          d_this->d_Sintra[i][g] = sinc(d_this->d_intraC2Wdist[i] * (Qmin + g * dQ));
        }
      }
    }
  }

  void NS::printSintra(ostream &os) {
    for(unsigned int i = 0; i < d_this->d_Sintra.size(); ++i) {
      for(unsigned int g = 0; g < d_this->d_Sintra[i].size(); ++g) {
        os << d_this->d_Sintra[i][g] << endl;
      }
      os << endl;
    }
  }

  void NS::printSinter(ostream &os) {
    for(unsigned int i = 0; i < d_this->d_Sinter.size(); ++i) {
      for(unsigned int g = 0; g < d_this->d_Sinter[i].size(); ++g) {
        os << d_this->d_Sinter[i][g] << endl;
      }
      os << endl;
    }
  }

  void NS::printS(std::ostream& os) {
    double Qmin = 2 * physConst.get_pi() * double(d_this->d_grid) /
            ((double(d_this->d_grid) - 0.5) * d_this->d_cut);
    double dQ = (d_this->d_Qmax - Qmin) / double(d_this->d_grid - 1);
    multimap<int, int>::iterator it;
    os << "#" << setw(19) << "Q";
    for(it = d_this->d_comb.begin(); it != d_this->d_comb.end(); ++it) {
      os << setw(20) << it->first+1 << "-" << it->second+1;
    }
    os << endl;
    int count = 0;
    for (int g = 0; g < d_this->d_grid; ++g) {
      os << setw(20) << Qmin + g * dQ;
      for (it = d_this->d_comb.begin(); it != d_this->d_comb.end(); ++it, ++count) {
        os << setw(20) << d_this->d_Sintra[count][g];
      }
      os << endl;
    }
  }

  void NS::calcSinter(void) {
    double Qmin = 2 * physConst.get_pi() * double(d_this->d_grid) /
            ((double(d_this->d_grid) - 0.5) * d_this->d_cut);
    double dQ = (d_this->d_Qmax - Qmin) / double(d_this->d_grid - 1);
    double dr = d_this->d_cut / d_this->d_grid;
    // loop over different RDFs
    for (unsigned int i = 0; i < d_this->d_rdf.size(); ++i) {
      // loop over the Q values
      for (int q = 0; q < d_this->d_grid; ++q) {
        double Q = Qmin + dQ * q;
        // loop over the centre to with distances of the RDF
        for (int r = 0; r < d_this->d_grid; ++r) {
          double R = (double(r) + 0.5) * dr;
          d_this->d_Sinter[i][q] += 4 * physConst.get_pi() * R * R * (d_this->d_rdf[i].rdf(r))
                  * sinc(R * Q) * dr;
        } // end of loop over RDF distances
      } // end of loop over Q valies
    } // end of loop over different RDFs
  }

  void NS::calcIntensity(void) {
    // make sure the intensity is empty
    for(int g = 0; g < d_this->d_grid; ++g) {
      d_this->d_intensity[g] = 0.0;
    }
    multimap<int, int>::const_iterator it;
    int count = 0;
    for(it = d_this->d_comb.begin(); it != d_this->d_comb.end(); ++it, ++count) {
      double bi = d_this->d_scattLen.find(it->first)->second;
      double bj = d_this->d_scattLen.find(it->second)->second;
      double bibj = bi * bj;
      double w_intra = d_this->d_weightIntra[count];
      double w_inter = d_this->d_weightInter[count];
      for(int g = 0; g < d_this->d_grid; ++g) {
        d_this->d_intensity[g] += w_intra * bibj * d_this->d_Sintra[count][g];
        d_this->d_intensity[g] += w_inter * bibj * d_this->d_Sinter[count][g];
      }
    }
  }

}
