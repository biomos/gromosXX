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

// gio_OutG96.cc
#include "OutG96.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <string>

#include "../gromos/Exception.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/Solvent.h"
#include "../gmath/Vec.h"
#include "../gmath/Matrix.h"
#include "../gcore/Box.h"
#include "../gcore/Remd.h"
#include "../args/Arguments.h"
#include "../utils/groTime.h"
#include "../utils/AtomSpecifier.h"
#include "OutCoordinates.h"

using gio::OutG96;
using namespace gcore;
using namespace std;
using namespace utils;

class gio::OutG96_i{
  friend class gio::OutG96;
  ostream &d_os;
  int d_count;
  int d_switch;
  OutG96_i(ostream &os):
    d_os(os), d_count(0)
  {d_switch = 0;}
  ~OutG96_i(){}

  void select(const string &thing);
  void writeTrajM(const Molecule &mol);
  void writeTrajV(const gcore::System &sys);
  void writeTrajS(const Solvent &sol);
  void writeBox(const Box &box);
  void writeTriclinicBox(const Box &box);
  void writeGenBox(const Box &box);
  void writeRemd(const Remd &remd);
  void writeAtomSpecifier(const AtomSpecifier & atoms);
};

OutG96::OutG96(ostream &os):
  OutCoordinates(),
  d_this(new OutG96_i(os)){}
OutG96::OutG96():
  OutCoordinates()
{d_this=0;}

OutG96::~OutG96(){
  if(d_this)delete d_this;
}

void OutG96::writeGenBox(const gcore::Box &box){
  d_this->writeGenBox(box);
}

void OutG96::writeTriclinicBox(const gcore::Box &box){
  d_this->writeTriclinicBox(box);
}

void OutG96::writeTitle(const string &title){
  d_this->d_os << "TITLE\n" << title << "\nEND\n";
}

void OutG96::writeTimestep(const int step, const double time)
{
  d_this->d_os.precision(9);
  d_this->d_os.setf(std::ios::fixed, std::ios::floatfield);
  
  d_this->d_os << "TIMESTEP\n"
	       << std::setw(15)
	       << step
	       << std::setw(20)
	       << time
               << "\n#if @time flag is used the value for step refers to the"
               << "\n#step-th configuration in the original trajectory file"
	       << "\nEND\n";
}

void OutG96::select(const string &thing){
  if (thing == "ALL"){
    d_this->d_switch = 1;
  } else if (thing =="SOLVENT"){
    d_this->d_switch = 2;
  } else if (thing == "SOLUTEV") {
    d_this->d_switch = 3;
  } else if (thing == "ALLV") {
    d_this->d_switch = 4;
  } else if (thing == "SOLVENTV") {
    d_this->d_switch = 5;
  } else {
    d_this->d_switch = 0;
  }
}

void OutG96::open(ostream &os){
  if(d_this){
    delete d_this;
  }
  d_this=new OutG96_i(os);
}

void OutG96::close(){
  if(d_this)delete d_this;
  d_this=0;
}

OutG96 &OutG96::operator<<(const gcore::System &sys){
  d_this->d_count=0;
  if(sys.hasRemd){
    d_this->d_os << "REMD\n";
    d_this->writeRemd(sys.remd());
    d_this->d_os << "END\n";
  }
  d_this->d_os << "POSITIONRED\n";
  if (d_this->d_switch == 0 || d_this->d_switch == 1 || 
      d_this->d_switch == 3 || d_this->d_switch == 4) {
    for(int i=0;i<sys.numMolecules();++i){
      d_this->writeTrajM(sys.mol(i));}
  }
  if (d_this->d_switch == 3 || d_this->d_switch == 4 || 
      d_this->d_switch == 5) {
    d_this->writeTrajV(sys);}
  if (d_this->d_switch == 1 || d_this->d_switch == 2 ||
      d_this->d_switch == 4 || d_this->d_switch == 5) {
    for(int i=0;i<sys.numSolvents();++i){
      d_this->writeTrajS(sys.sol(i));}
  }
  d_this->d_os << "END\n";

  if (args::Arguments::outG96) {
    switch (sys.box().boxformat()) {
      case gcore::Box::box96 :
        d_this->d_os << "BOX\n";
        d_this->writeBox(sys.box());
        d_this->d_os << "END\n";
        break;
      case gcore::Box::triclinicbox :
        d_this->d_os << "TRICLINICBOX\n";
        d_this->writeTriclinicBox(sys.box());
        d_this->d_os << "END\n";
        break;
      case gcore::Box::genbox :
        d_this->d_os << "GENBOX\n";
        d_this->writeGenBox(sys.box());
        d_this->d_os << "END\n";
        break;
      default:
        throw gromos::Exception("OutG96", "Don't know how to handle boxformat");
    }
  } else {
    // in GXX there is only one single format called GENBOX
    d_this->d_os << "GENBOX\n";
    d_this->writeGenBox(sys.box());
    d_this->d_os << "END\n";      
  }
  
  return *this;
}

OutG96 &OutG96::operator<<(const utils::AtomSpecifier & atoms){
  const System & sys = *(atoms.sys());
  d_this->d_count=0;
  if(sys.hasRemd){
    d_this->d_os << "REMD\n";
    d_this->writeRemd(sys.remd());
    d_this->d_os << "END\n";
  }
  d_this->writeAtomSpecifier(atoms);

  if (args::Arguments::outG96) {
    switch (sys.box().boxformat()) {
      case gcore::Box::box96 :
        d_this->d_os << "BOX\n";
        d_this->writeBox(sys.box());
        d_this->d_os << "END\n";
        break;
      case gcore::Box::triclinicbox :
        d_this->d_os << "TRICLINICBOX\n";
        d_this->writeTriclinicBox(sys.box());
        d_this->d_os << "END\n";
        break;
      case gcore::Box::genbox :
        d_this->d_os << "GENBOX\n";
        d_this->writeGenBox(sys.box());
        d_this->d_os << "END\n";
        break;
      default:
        throw gromos::Exception("OutG96", "Don't know how to handle boxformat");
    }
  } else {
    // in GXX there is only one single format called GENBOX
    d_this->d_os << "GENBOX\n";
    d_this->writeGenBox(sys.box());
    d_this->d_os << "END\n";
  }

  return *this;
}

void gio::OutG96_i::writeRemd(const Remd &remd)
{
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.precision(1);
  d_os << setw(15) << remd.id() 
       << setw(10) << remd.run()
       << setprecision(1) << setw(10) << remd.temperature()
       << setprecision(6) << setw(10) << remd.lambda()
       << "\n"
       << setw(15) << remd.Ti()
       << setw(10) << remd.li()
       << setw(10) << remd.Tj()
       << setw(10) << remd.lj()
       << setw(10) << remd.reeval()
       << "\n";
}

void gio::OutG96_i::writeTrajM(const Molecule &mol){
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.precision(9);
  for (int i=0;i<mol.numPos();++i){
    ++d_count;
    d_os << setw(15) << mol.pos(i)[0]
	 << setw(15) << mol.pos(i)[1]
	 << setw(15) << mol.pos(i)[2]<< endl;
    if(!(d_count%10))
      d_os << "#" << setw(10)<<d_count<<endl;
  } 
}
void gio::OutG96_i::writeTrajV(const gcore::System &sys){
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.precision(9);
  for (int i=0;i<sys.vas().numVirtualAtoms();++i){
    ++d_count;
    d_os << setw(15) << sys.vas().atom(i).pos()[0]
	 << setw(15) << sys.vas().atom(i).pos()[1]
	 << setw(15) << sys.vas().atom(i).pos()[2]<< endl;
    if(!(d_count%10))
      d_os << "#" << setw(10)<<d_count<<endl;
  } 
}
void gio::OutG96_i::writeTrajS(const Solvent &sol){
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.precision(9);
  for (int i=0;i<sol.numPos();++i){
    ++d_count;
    d_os << setw(15) << sol.pos(i)[0]
	 << setw(15) << sol.pos(i)[1]
	 << setw(15) << sol.pos(i)[2]<< endl;
    if(!(d_count%10))
      d_os << "#" << setw(10)<<d_count<<endl;
  } 
}


void gio::OutG96_i::writeBox(const Box &box){
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.precision(9);

  d_os << setw(15) << box.K()[0]
       << setw(15) << box.L()[1]
       << setw(15) << box.M()[2] << endl;
}

void gio::OutG96_i::writeTriclinicBox(const Box &box){
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.precision(9);

  d_os << setw(8) << box.ntb() << endl;
  for(int i=0; i<3; ++i){
    d_os << setw(15) << box.K()[i] 
	 << setw(15) << box.L()[i]
	 << setw(15) << box.M()[i] << endl;
  }
}

void gio::OutG96_i::writeGenBox(const Box &box){
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.precision(9);
  const double k=box.K().abs();
  const double l=box.L().abs();
  const double m=box.M().abs();
  d_os << setw(8) << box.ntb() << endl;
  if(box.ntb()==gcore::Box::vacuum)
    d_os << setw(15) << 0.0 << setw(15) << 0.0 << setw(15) << 0.0 << endl
	 << setw(15) << 0.0 << setw(15) << 0.0 << setw(15) << 0.0 << endl
          << setw(15) << 0.0 << setw(15) << 0.0 << setw(15) << 0.0 << endl
          << setw(15) << box.X() << setw(15) << box.Y() << setw(15) << box.Z() << endl;
  else{
    d_os << setw(15) << k
	 << setw(15) << l
	 << setw(15) << m << endl;
    d_os << setw(15) << acos(box.L().dot(box.M())/(l*m))*180/M_PI
	 << setw(15) << acos(box.K().dot(box.M())/(k*m))*180/M_PI
	 << setw(15) << acos(box.K().dot(box.L())/(k*l))*180/M_PI << endl;

    // calculate the Euler rotation angles as described in Phils manuscript:
    // "GROMOS01: Description of the changes", Philippe Huenenberger, October 5, 2004
    gmath::Vec x = box.K().normalize();
    gmath::Vec y = (box.L() - (box.L().dot(x) * x)).normalize();
    gmath::Vec z = x.cross(y);

    gmath::Matrix R_(x, y, z);
    double R11R21 = R_(0,0) * R_(0,0) + R_(1,0) * R_(1,0);
    double theta, psi, phi;
    if(R11R21 == 0.0) {
        int sign = 1;
        if(R_(2,0)<0) sign = -1;
        theta = -sign*M_PI/2;
        psi = 0.0;
        sign = 1;
        if(R_(0,1)<0) sign = -1;
        phi = -sign*acos(R_(1,1));
    } else {
        int sign =1;
        if(R_(2,0)<0) sign = -1;
        theta = -sign*acos(sqrt(R_(0,0)*R_(0,0)+R_(1,0)*R_(1,0)));
        sign = 1;
        if((R_(2,1)/cos(theta))<0) sign = -1;
        psi = sign*acos(R_(2,2)/cos(theta));
        sign = 1;
        if((R_(1,0)/cos(theta))<0) sign = -1;
        phi = sign*acos(R_(0,0)/cos(theta));
    }
    
    d_os << setw(15) << phi/M_PI*180
	 << setw(15) << theta/M_PI*180
	 << setw(15) << psi/M_PI*180 << endl;

    d_os << setw(15) << box.X()
            << setw(15) << box.Y()
            << setw(15) << box.Z() << endl;
  }
}

void gio::OutG96_i::writeAtomSpecifier(const AtomSpecifier& atoms) {
  d_os << "POSITIONRED" << endl;
  d_os.setf(ios::fixed, ios::floatfield);
  d_os.precision(9);
  d_os << "# selected " << atoms.size() << " atoms" << endl;
  for (unsigned int i = 0; i < atoms.size(); ++i) {
    d_os << setw(15) << atoms.pos(i)[0]
            << setw(15) << atoms.pos(i)[1]
            << setw(15) << atoms.pos(i)[2]
            << "  # ";
    if (atoms.mol(i) < 0)
      d_os << "s";
    else 
      d_os << atoms.mol(i) + 1;
    d_os << ":" << atoms.atom(i) + 1 << endl;
    if (!((i + 1) % 10))
      d_os << "#" << setw(10) << i + 1 << endl;
  }
  d_os << "END" << endl;
}
