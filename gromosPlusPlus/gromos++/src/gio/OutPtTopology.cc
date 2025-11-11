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
#include "OutPtTopology.h"

#include <cstddef>
#include <utility>
#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <iomanip>
#include <cassert>

#include "OutTopology.h"
#include "../gcore/Bond.h"
#include "../gcore/Angle.h"
#include "../gcore/Dihedral.h"
#include "../gcore/CrossDihedral.h"
#include "../gcore/Improper.h"
#include "../gcore/AtomPair.h"
#include "../gcore/System.h"
#include "../gcore/PtTopology.h"
#include "../gcore/System.h"
#include "../gcore/LinearTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gromos/Exception.h"

using namespace std;
using namespace gcore;

namespace gio{
  OutPtTopology::OutPtTopology(ostream &os) :
  d_title("a perturbation topology"), d_os(os) {
  }
 
  OutPtTopology::~OutPtTopology() {}
  
  void OutPtTopology::setTitle(const string &title) {
    d_title = title;
  }
  
  void OutPtTopology::write(const PtTopology & pt, System * sys, int a, int b) {
    LinearTopology * top = NULL;
    if (sys != NULL) {
      top = new LinearTopology(*sys);
    }
    
    if (pt.numPt() != 2) {
      throw gromos::Exception("OutPtTopology", "Can only write a topology with "
              "two states.");
    }
    
    // write the title
    d_os << "TITLE" << endl
         << d_title << endl
         << "END" << endl;
    
    // atoms
    d_os << "PERTATOMPARAM" << endl
         << "# number of perturbed atoms" << endl
         << setw(5) << pt.numAtoms() << endl
         << "#  NR  RES  NAME IAC(A)  MASS(A) CHARGE(A) IAC(B)  MASS(B) "
            "CHARGE(B)       ALJ       ACRF" << endl;
    for (int i = 0; i < pt.numAtoms(); ++i) {
      d_os.precision(5);
      d_os.setf(ios::fixed, ios::floatfield);
      d_os << setw(5) << pt.atomNum(i)+1 << setw(5) << pt.residueNum(i)+1 << " "
              << setw(5) << pt.atomName(i) << setw(4) << pt.iac(i,a)+1 << setw(11)
              << pt.mass(i,a) << setw(11) << pt.charge(i,a)
              << setw(4) << pt.iac(i,b)+1 << setw(11) << pt.mass(i,b) << setw(11) << pt.charge(i,b);
      d_os.precision(2);
      d_os << setw(11) << pt.alphaLJ(i) << setw(11) << pt.alphaCRF(i) << endl;
    }
    d_os << "END" << endl;

    // polarisation
    if (pt.hasPolarisationParameters()) {
      d_os << "PERTATOMPARAM" << endl
              << "# number of perturbed atoms" << endl
              << setw(5) << pt.numAtoms() << endl
              << "#  NR  RES  NAME IAC(A)  MASS(A) CHARGE(A) IAC(B)  MASS(B) "
              "CHARGE(B)       ALJ       ACRF" << endl;
      for (int i = 0; i < pt.numAtoms(); ++i) {
        d_os.precision(5);
        d_os.setf(ios::fixed, ios::floatfield);
        d_os << setw(5) << pt.atomNum(i)+1 << setw(5) << pt.residueNum(i)+1 << " "
             << setw(5) << pt.atomName(i)
             << setw(11) << pt.polarisability(i,a) << setw(11) << pt.dampingLevel(i,a)
             << setw(11) << pt.polarisability(i,b) << setw(11) << pt.dampingLevel(i,b) << endl;
      }
      d_os << "END" << endl;
    }
    
    // atom pairs
    if (pt.atompairs(a).size() != pt.atompairs(b).size())
      throw gromos::Exception("OutPtTopology", "Number of atom pairs are "
              "different in the two states.");
    
    d_os << "PERTATOMPAIR" << endl
         << "# number of perturbed atom pairs" << endl
         << pt.atompairs(a).size() << endl
         << "#    i     j  i(A)  i(B)" << endl;
    {
      set<AtomPairParam>::const_iterator it_a = pt.atompairs(a).begin(),
              it_b = pt.atompairs(b).begin(), to = pt.atompairs(a).end();
      for(; it_a != to; ++it_a, ++it_b) {
        d_os << setw(6) << (*it_a)[0]+1 << setw(6) << (*it_a)[1]+1
             << setw(6) << it_a->param() << setw(6) << it_b->param() << endl;
      }
    }
    d_os << "END" << endl;

    // bonds
    if (pt.bonds(a).size() != pt.bonds(b).size())
      throw gromos::Exception("OutPtTopology", "Number of bonds are "
              "different in the two states.");
    
    set<pair<Bond, Bond> > bonds;
    set<pair<Bond, Bond> > bondsH;
    
    if (top != NULL) { // partition into hydrogens and non hydrogens
      set<Bond>::const_iterator it_a = pt.bonds(a).begin(),
              to = pt.bonds(a).end(),
              it_b = pt.bonds(b).begin();
      for (; it_a != to; ++it_a, ++it_b) {
        if (top->atoms()[(*it_a)[0]].mass() == 1.008 || 
            top->atoms()[(*it_a)[1]].mass() == 1.008) {
          bondsH.insert(pair<Bond, Bond>(*it_a, *it_b));
        } else {
          bonds.insert(pair<Bond, Bond>(*it_a, *it_b));
        }
      }
    } else { // don't partition into H, non-H
      set<Bond>::const_iterator it_a = pt.bonds(a).begin(),
              to = pt.bonds(a).end(),
              it_b = pt.bonds(b).begin();
      for (; it_a != to; ++it_a, ++it_b) {
        bonds.insert(pair<Bond, Bond>(*it_a, *it_b));
      }      
    } // end partitioning
    
    d_os << "PERTBONDSTRETCH" << endl
         << "# number of perturbed bonds" << endl
         << bonds.size() << endl
         << "#    i     j t(A) t(B)" << endl;
    for(set<pair<Bond,Bond> >::iterator it = bonds.begin(), to = bonds.end();
            it != to; ++it) {
      d_os << setw(6) << it->first[0]+1 << setw(6) << it->first[1]+1
           << setw(5) << it->first.type()+1 << setw(5) << it->second.type()+1 << endl;
    }
    d_os << "END" << endl;
    
    d_os << "PERTBONDSTRETCHH" << endl
         << "# number of perturbed bonds" << endl
         << bondsH.size() << endl
         << "#    i     j t(A) t(B)" << endl;
    for(set<pair<Bond,Bond> >::iterator it = bondsH.begin(), to = bondsH.end();
            it != to; ++it) {
      d_os << setw(6) << it->first[0]+1 << setw(6) << it->first[1]+1
           << setw(5) << it->first.type()+1 << setw(5) << it->second.type()+1 << endl;
    }
    d_os << "END" << endl;
    
    // angles
    if (pt.angles(a).size() != pt.angles(b).size())
      throw gromos::Exception("OutPtTopology", "Number of bond angles are "
              "different in the two states.");
    set<pair<Angle, Angle> > angles;
    set<pair<Angle, Angle> > anglesH;

    if (top != NULL) { // partition into hydrogens and non hydrogens
      set<Angle>::const_iterator it_a = pt.angles(a).begin(),
              to = pt.angles(a).end(),
              it_b = pt.angles(b).begin();
      for (; it_a != to; ++it_a, ++it_b) {
        if (top->atoms()[(*it_a)[0]].mass() == 1.008 || 
            top->atoms()[(*it_a)[1]].mass() == 1.008 ||
            top->atoms()[(*it_a)[2]].mass() == 1.008) {
          anglesH.insert(pair<Angle, Angle>(*it_a, *it_b));
        } else {
          angles.insert(pair<Angle, Angle>(*it_a, *it_b));
        }
      }
    } else { // don't partition into H, non-H
      set<Angle>::const_iterator it_a = pt.angles(a).begin(),
              to = pt.angles(a).end(),
              it_b = pt.angles(b).begin();
      for (; it_a != to; ++it_a, ++it_b) {
        angles.insert(pair<Angle, Angle>(*it_a, *it_b));
      }      
    } // end partitioning
    
    d_os << "PERTBONDANGLE" << endl
         << "# number of perturbed bond angles" << endl
         << angles.size() << endl
         << "#    i     j     k t(A) t(B)" << endl;
    for(set<pair<Angle,Angle> >::iterator it = angles.begin(), to = angles.end();
            it != to; ++it) {
      d_os << setw(6) << it->first[0]+1 << setw(6) << it->first[1]+1 
           << setw(6) << it->first[2]+1
           << setw(5) << it->first.type()+1 << setw(5) << it->second.type()+1 << endl;
    }
    d_os << "END" << endl;
    
    d_os << "PERTBONDANGLEH" << endl
         << "# number of perturbed bond angles" << endl
         << anglesH.size() << endl
         << "#    i     j     k t(A) t(B)" << endl;
    for(set<pair<Angle,Angle> >::iterator it = anglesH.begin(), to = anglesH.end();
            it != to; ++it) {
      d_os << setw(6) << it->first[0]+1 << setw(6) << it->first[1]+1 
           << setw(6) << it->first[2]+1
           << setw(5) << it->first.type()+1 << setw(5) << it->second.type()+1 << endl;
    }
    d_os << "END" << endl;
    
    // impropers
    if (pt.impropers(a).size() != pt.impropers(b).size())
      throw gromos::Exception("OutPtTopology", "Number of improper dihedrals are "
              "different in the two states.");
    set<pair<Improper, Improper> > impropers;
    set<pair<Improper, Improper> > impropersH;

    if (top != NULL) { // partition into hydrogens and non hydrogens
      set<Improper>::const_iterator it_a = pt.impropers(a).begin(),
              to = pt.impropers(a).end(),
              it_b = pt.impropers(b).begin();
      for (; it_a != to; ++it_a, ++it_b) {
        if (top->atoms()[(*it_a)[0]].mass() == 1.008 || 
            top->atoms()[(*it_a)[1]].mass() == 1.008 ||
            top->atoms()[(*it_a)[2]].mass() == 1.008 ||
            top->atoms()[(*it_a)[3]].mass() == 1.008) {
          impropersH.insert(pair<Improper, Improper>(*it_a, *it_b));
        } else {
          impropers.insert(pair<Improper, Improper>(*it_a, *it_b));
        }
      }
    } else { // don't partition into H, non-H
      set<Improper>::const_iterator it_a = pt.impropers(a).begin(),
              to = pt.impropers(a).end(),
              it_b = pt.impropers(b).begin();
      for (; it_a != to; ++it_a, ++it_b) {
        impropers.insert(pair<Improper, Improper>(*it_a, *it_b));
      }      
    } // end partitioning
    
    d_os << "PERTIMPROPERDIH" << endl
         << "# number of perturbed improper dihedrals" << endl
         << impropers.size() << endl
         << "#    i     j     k     l t(A) t(B)" << endl;
    for(set<pair<Improper,Improper> >::iterator it = impropers.begin(), to = impropers.end();
            it != to; ++it) {
      d_os << setw(6) << it->first[0]+1 << setw(6) << it->first[1]+1 
           << setw(6) << it->first[2]+1 << setw(6) << it->first[3]+1
           << setw(5) << it->first.type()+1 << setw(5) << it->second.type()+1 << endl;
    }
    d_os << "END" << endl;
    
    d_os << "PERTIMPROPERDIHH" << endl
         << "# number of perturbed improper dihedrals" << endl
         << impropersH.size() << endl
         << "#    i     j     k     l t(A) t(B)" << endl;
    for(set<pair<Improper,Improper> >::iterator it = impropersH.begin(), to = impropersH.end();
            it != to; ++it) {
      d_os << setw(6) << it->first[0]+1 << setw(6) << it->first[1]+1 
           << setw(6) << it->first[2]+1 << setw(6) << it->first[3]+1
           << setw(5) << it->first.type()+1 << setw(5) << it->second.type()+1 << endl;
    }
    d_os << "END" << endl;
    
    // dihedral
    if (pt.dihedrals(a).size() != pt.dihedrals(b).size())
      throw gromos::Exception("OutPtTopology", "Number of dihedrals are "
              "different in the two states.");
    set<pair<Dihedral, Dihedral> > dihedrals;
    set<pair<Dihedral, Dihedral> > dihedralsH;
    
    if (top != NULL) { // partition into hydrogens and non hydrogens
      set<Dihedral>::const_iterator it_a = pt.dihedrals(a).begin(),
              to = pt.dihedrals(a).end(),
              it_b = pt.dihedrals(b).begin();
      for (; it_a != to; ++it_a, ++it_b) {
        if (top->atoms()[(*it_a)[0]].mass() == 1.008 || 
            top->atoms()[(*it_a)[1]].mass() == 1.008 ||
            top->atoms()[(*it_a)[2]].mass() == 1.008 ||
            top->atoms()[(*it_a)[3]].mass() == 1.008) {
          dihedralsH.insert(pair<Dihedral, Dihedral>(*it_a, *it_b));
        } else {
          dihedrals.insert(pair<Dihedral, Dihedral>(*it_a, *it_b));
        }
      }
    } else { // don't partition into H, non-H
      set<Dihedral>::const_iterator it_a = pt.dihedrals(a).begin(),
              to = pt.dihedrals(a).end(),
              it_b = pt.dihedrals(b).begin();
      for (; it_a != to; ++it_a, ++it_b) {
        dihedrals.insert(pair<Dihedral, Dihedral>(*it_a, *it_b));
      }      
    } // end partitioning
    
    d_os << "PERTPROPERDIH" << endl
         << "# number of perturbed dihedrals" << endl
         << dihedrals.size() << endl
         << "#    i     j     k     l t(A) t(B)" << endl;
    for(set<pair<Dihedral,Dihedral> >::iterator it = dihedrals.begin(), to = dihedrals.end();
            it != to; ++it) {
      d_os << setw(6) << it->first[0]+1 << setw(6) << it->first[1]+1 
           << setw(6) << it->first[2]+1 << setw(6) << it->first[3]+1
           << setw(5) << it->first.type()+1 << setw(5) << it->second.type()+1 << endl;
    }
    d_os << "END" << endl;
    
    d_os << "PERTPROPERDIHH" << endl
         << "# number of perturbed dihedrals" << endl
         << dihedralsH.size() << endl
         << "#    i     j     k     l t(A) t(B)" << endl;
    for(set<pair<Dihedral,Dihedral> >::iterator it = dihedralsH.begin(), to = dihedralsH.end();
            it != to; ++it) {
      d_os << setw(6) << it->first[0]+1 << setw(6) << it->first[1]+1 
           << setw(6) << it->first[2]+1 << setw(6) << it->first[3]+1
           << setw(5) << it->first.type()+1 << setw(5) << it->second.type()+1 << endl;
    }
    d_os << "END" << endl;

    // crossdihedral
    if (pt.crossdihedrals(a).size() != pt.crossdihedrals(b).size())
      throw gromos::Exception("OutPtTopology", "Number of cross dihedrals are "
              "different in the two states.");
    set<pair<CrossDihedral, CrossDihedral> > crossdihedrals;
    set<pair<CrossDihedral, CrossDihedral> > crossdihedralsH;

    if (top != NULL) { // partition into hydrogens and non hydrogens
      set<CrossDihedral>::const_iterator it_a = pt.crossdihedrals(a).begin(),
              to = pt.crossdihedrals(a).end(),
              it_b = pt.crossdihedrals(b).begin();
      for (; it_a != to; ++it_a, ++it_b) {
        if (top->atoms()[(*it_a)[0]].mass() == 1.008 ||
            top->atoms()[(*it_a)[1]].mass() == 1.008 ||
            top->atoms()[(*it_a)[2]].mass() == 1.008 ||
            top->atoms()[(*it_a)[3]].mass() == 1.008 ||
            top->atoms()[(*it_a)[4]].mass() == 1.008 ||
            top->atoms()[(*it_a)[5]].mass() == 1.008 ||
            top->atoms()[(*it_a)[6]].mass() == 1.008 ||
            top->atoms()[(*it_a)[7]].mass() == 1.008) {
          crossdihedralsH.insert(pair<CrossDihedral, CrossDihedral>(*it_a, *it_b));
        } else {
          crossdihedrals.insert(pair<CrossDihedral, CrossDihedral>(*it_a, *it_b));
        }
      }
    } else { // don't partition into H, non-H
      set<CrossDihedral>::const_iterator it_a = pt.crossdihedrals(a).begin(),
              to = pt.crossdihedrals(a).end(),
              it_b = pt.crossdihedrals(b).begin();
      for (; it_a != to; ++it_a, ++it_b) {
        crossdihedrals.insert(pair<CrossDihedral, CrossDihedral>(*it_a, *it_b));
      }
    } // end partitioning

    // Currently MD++ cannot perturb cross dihedrals. 
    // So write this block only when it is really necessary
    if(crossdihedrals.size()){
	
      d_os << "PERTCROSSDIH" << endl
	   << "# WARNING: MD++ may not know how to handle perturbed cross dihedrals"
	   << "# number of perturbed cross dihedrals" << endl
	   << crossdihedrals.size() << endl
	   << "#    a     b     c     d     e     f     g     h t(A) t(B)" << endl;
      for(set<pair<CrossDihedral,CrossDihedral> >::iterator it = crossdihedrals.begin(), to = crossdihedrals.end();
	  it != to; ++it) {
	d_os << setw(6) << it->first[0]+1 << setw(6) << it->first[1]+1
	     << setw(6) << it->first[2]+1 << setw(6) << it->first[3]+1
	     << setw(6) << it->first[4]+1 << setw(6) << it->first[5]+1
	     << setw(6) << it->first[6]+1 << setw(6) << it->first[7]+1
	     << setw(5) << it->first.type()+1 << setw(5) << it->second.type()+1 << endl;
      }
      d_os << "END" << endl;
    }
    if(crossdihedralsH.size()){
	
      d_os << "PERTCROSSDIH" << endl
	   << "# WARNING: MD++ may not know how to handle perturbed cross dihedrals"
	   << "# number of perturbed cross dihedrals" << endl
	   << crossdihedralsH.size() << endl
	   << "#    a     b     c     d     e     f     g     h t(A) t(B)" << endl;
      for(set<pair<CrossDihedral,CrossDihedral> >::iterator it = crossdihedralsH.begin(), to = crossdihedralsH.end();
	  it != to; ++it) {
	d_os << setw(6) << it->first[0]+1 << setw(6) << it->first[1]+1
	     << setw(6) << it->first[2]+1 << setw(6) << it->first[3]+1
	     << setw(6) << it->first[4]+1 << setw(6) << it->first[5]+1
	     << setw(6) << it->first[6]+1 << setw(6) << it->first[7]+1
	     << setw(5) << it->first.type()+1 << setw(5) << it->second.type()+1 << endl;
      }
      d_os << "END" << endl;
    }
    
    if (top != NULL)
      delete top;
  }
  
  void OutPtTopology::write_multiple(const PtTopology & pt) {    
    // write the title
    d_os << "TITLE" << endl
         << d_title << endl
         << "END" << endl;
    
    d_os << "MPERTATOM" << endl
         << "# NJLA: number of perturbed atoms" << endl
         << "# NPTB: number of states" << endl
         << setw(5) << pt.numAtoms() << setw(5) << pt.numPt() << endl
         << "# identifiers of the states" << endl;
         
    // the name of the states
    for(int i = 0; i < pt.numPt(); ++i)
      d_os << setw(10) << pt.pertName(i);
    d_os << endl;
    d_os << "#  NR  NAME IAC(1) CHARGE(1) ...  IAC(n) CHARGE(n)  ALPHLJ  ALPHCRF" << endl;
    for(int i = 0; i < pt.numAtoms(); ++i) {
      d_os.precision(5);
      d_os.setf(ios::fixed, ios::floatfield);
      d_os << setw(5) << pt.atomNum(i)+1 << " " << setw(5) << pt.atomName(i);
      for(int j = 0; j < pt.numPt(); ++j) {
        d_os << setw(4) << pt.iac(i,j)+1 << setw(11) << pt.charge(i,j);
      }
      d_os << setw(11) << pt.alphaLJ(i) << setw(11) << pt.alphaCRF(i) << endl;
    } // for atoms
    d_os << "END" << endl;
    
    for(int j = 0; j < pt.numPt(); ++j) {
      if (pt.bonds(j).size() || pt.angles(j).size() || pt.impropers(j).size() ||
              pt.atompairs(j).size() || pt.dihedrals(j).size() || pt.crossdihedrals(j).size()) {
        throw gromos::Exception("OutPtTopology", "There are bonded parameters "
                "or exclusion changes in a multiple perturbation topology.");
      }
    } // for states
  }
} // namespace
