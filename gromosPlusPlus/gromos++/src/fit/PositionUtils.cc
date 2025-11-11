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

// fit_PositionUtils.cc
#include "PositionUtils.h"

#include <math.h>

#include "Reference.h"
#include "../gmath/Vec.h"
#include "../gmath/Matrix.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Solvent.h"
#include "../gcore/AtomTopology.h"
#include "../utils/AtomSpecifier.h"


using namespace gmath;
using namespace gcore;
using fit::PositionUtils;

Vec PositionUtils::com(const gcore::System &sys){
  // calculate the center of mass of the molecule
  double totalMass=0;
  Vec cm;

  for(int m=0;m<sys.numMolecules();++m)
    for (int i=0;i < sys.mol(m).numAtoms(); i++) {
      cm += sys.mol(m).pos(i) * sys.mol(m).topology().atom(i).mass();
      totalMass += sys.mol(m).topology().atom(i).mass();
    }

  cm = (1.0/totalMass)*cm;
  
  return cm;
}

Vec PositionUtils::cog(const gcore::System &sys){
  
  // calculate the center of mass of the molecule
  Vec cm;
  int atoms=0;

  for(int m=0;m<sys.numMolecules();++m)
    for (int i=0;i < sys.mol(m).numAtoms(); i++) {
      cm = cm + sys.mol(m).pos(i);
      ++atoms;
    }

  cm = (1.0/double(atoms))*cm;
  
  return cm;
}

Vec PositionUtils::com(const System &sys, const Reference &ref){
  double totalMass=0;
  Vec cm;

  for(int m=0;m<sys.numMolecules();++m)
    for (int i=0;i < sys.mol(m).numAtoms(); i++) {
      cm += 
	ref.weight(m,i)
	* sys.mol(m).topology().atom(i).mass()
	* sys.mol(m).pos(i) ;

      totalMass += ref.weight(m,i) *
	sys.mol(m).topology().atom(i).mass();

    }

  cm = (1.0/totalMass)*cm;
  return cm;
}

Vec PositionUtils::com(const System &sys, utils::AtomSpecifier & atoms){
  double totalMass=0;
  Vec cm;

  for(unsigned int a=0; a<atoms.size(); ++a){
    cm += 
      atoms.mass(a) * atoms.pos(a) ;
    
    totalMass += atoms.mass(a);
    
  }

  cm /= totalMass;
  return cm;
}

Vec PositionUtils::com_v(const System &sys, utils::AtomSpecifier & atoms){
  double totalMass=0;
  Vec cm_v;

  for(unsigned int a=0; a<atoms.size(); ++a){
    cm_v += 
      atoms.mass(a) * atoms.vel(a) ;
    
    totalMass += atoms.mass(a);
    
  }

  cm_v /= totalMass;
  return cm_v;
}

Vec PositionUtils::cog(const System &sys, const Reference &ref){
  Vec cg;

  for(int m=0;m<sys.numMolecules();++m)
    for (int i=0;i < sys.mol(m).numAtoms(); i++) {

     cg += 
	ref.weight(m,i)
	* sys.mol(m).pos(i) ;
    }

  return cg;
}

Vec PositionUtils::cog(const System &sys, utils::AtomSpecifier & atoms){
  Vec cg;

  for(unsigned int a=0;a<atoms.size();++a)
    cg += atoms.pos(a);
  
  cg /= atoms.size();
  
  return cg;
}

void PositionUtils::translate(gcore::System *sys, const gmath::Vec &v){
  for(int m=0;m<sys->numMolecules();++m)
    for(int i=0;i<sys->mol(m).numAtoms();++i)
      sys->mol(m).pos(i)=sys->mol(m).pos(i)+v;
  for(int j=0;j<sys->sol(0).numPos(); ++j)
    sys->sol(0).pos(j)=sys->sol(0).pos(j)+v;
}

void PositionUtils::translate(gcore::Molecule & mol, const gmath::Vec &v){
  for(int i=0; i < mol.numAtoms(); ++i)
    mol.pos(i) = mol.pos(i) + v;
}

void PositionUtils::rotate(gcore::System *sys, const gmath::Matrix &mat){
  for(int j=0;j<sys->numMolecules();++j)
    for(int i=0;i<sys->mol(j).numAtoms();++i)
      sys->mol(j).pos(i)=mat*sys->mol(j).pos(i);
  for(int j=0;j<sys->sol(0).numPos(); ++j)
    sys->sol(0).pos(j)= mat*sys->sol(0).pos(j);
  
}

void PositionUtils::rotate(gcore::Molecule & mol, const gmath::Matrix &mat){

  for(int i=0; i < mol.numAtoms(); ++i)
    mol.pos(i) = mat * mol.pos(i);
}

gmath::Matrix PositionUtils::rotateAround(gmath::Vec v, double a)
{
  // If speed is ever an issue for this function, we could of course
  // do the multiplications on paper and implement just the result.

  // calculate the angle in radians
  a*=M_PI/180.0;
  // first we create the matrix that brings v along z.
  // this is the product of two matrices
  // 1. around the x-axis by theta
  //    with sin(theta) = v[2]/r_yz; cos(theta) = align[2]/r_yz
  // 2. around the y-axis by phi
  //    with sin(phi) = align[0]/r; cos(phi) = r_yz / r

  double r=v.abs();
  double r_yz = sqrt(v[1]*v[1] + v[2]*v[2]);
  gmath::Matrix rot1(Vec( r_yz / r         ,  0         , v[0]/r ),
                     Vec(-v[0]*v[1]/r/r_yz ,  v[2]/r_yz , v[1]/r ),
                     Vec(-v[0]*v[2]/r/r_yz , -v[1]/r_yz , v[2]/r ));

  // we also need its transpose to rotate back
  gmath::Matrix rot2(3,3);
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      rot2(i,j)=rot1(j,i);
  
  // and a Matrix that rotates around z by the angle a
  double cosd=cos(a);
  double sind=sin(a);
  gmath::Matrix m3(Vec(cosd, sind, 0.0),
                   Vec(-sind, cosd, 0.0),
                   Vec(0.0, 0.0, 1.0));
  
  // The matrix that we want is now rot2 * m3 * rot1
  return rot2*m3*rot1;
}

Vec PositionUtils::shiftToCom(gcore::System *sys){
  Vec v= -com(*sys);
  translate(sys,v);
  return v;
}

Vec PositionUtils::shiftToCom(gcore::System *sys, const Reference &ref){
  Vec v=-com(*sys,ref);
  translate(sys, v);
  return v;
}

Vec PositionUtils::shiftToCom(gcore::System *sys, utils::AtomSpecifier & atoms){
  Vec v=-com(*sys, atoms);
  translate(sys, v);
  return v;
}

Vec PositionUtils::shiftToCog(gcore::System *sys){
  Vec v=-cog(*sys);
  translate(sys,v);
  return v;
}

Vec PositionUtils::shiftToCog(gcore::System *sys, const Reference &ref){
  Vec v= -cog(*sys, ref);
  translate(sys, v);
  return v;
}

Vec PositionUtils::shiftToCog(gcore::System *sys, utils::AtomSpecifier & atoms){
  Vec v= -cog(*sys, atoms);
  translate(sys, v);
  return v;
}

Vec PositionUtils::getmaxcoordinates(gcore::System *sys, bool heavy){
    
     double maxX = -100000000;
     double maxY = -100000000;
     double maxZ = -100000000;

     
      for(int i=0;i<sys->numMolecules();++i) {
       for(int j=0;j<sys->mol(i).numAtoms();++j) {
        Vec tmp = sys->mol(i).pos(j);
	if (heavy) {
	  if (sys->mol(i).topology().atom(j).mass() != 1.00800) {           
           maxX = ((tmp[0] > maxX) ? (tmp[0]) : (maxX));
           maxY = ((tmp[1] > maxY) ? (tmp[1]) : (maxY));
           maxZ = ((tmp[2] > maxZ) ? (tmp[2]) : (maxZ));
	  }
	}
        else {
         maxX = ((tmp[0] > maxX) ? (tmp[0]) : (maxX));
         maxY = ((tmp[1] > maxY) ? (tmp[1]) : (maxY));
         maxZ = ((tmp[2] > maxZ) ? (tmp[2]) : (maxZ));
	}
       }
      }

    
  return Vec (maxX,maxY,maxZ);
}
Vec PositionUtils::getmincoordinates(gcore::System *sys, bool heavy){
    
     double minX = 100000000;
     double minY = 100000000;
     double minZ = 100000000;

     
     for(int i=0;i<sys->numMolecules();++i) {
       for(int j=0;j<sys->mol(i).numAtoms();++j) {
       Vec tmp = sys->mol(i).pos(j);
	if (heavy) {
	  if (sys->mol(i).topology().atom(j).mass() != 1.00800) {           
           minX = ((tmp[0] < minX) ? (tmp[0]) : (minX));
           minY = ((tmp[1] < minY) ? (tmp[1]) : (minY));
           minZ = ((tmp[2] < minZ) ? (tmp[2]) : (minZ));
	  }
	}
        else {
         minX = ((tmp[0] < minX) ? (tmp[0]) : (minX));
         minY = ((tmp[1] < minY) ? (tmp[1]) : (minY));
         minZ = ((tmp[2] < minZ) ? (tmp[2]) : (minZ));
	}
       }
     }

  return Vec (minX,minY,minZ);
}
