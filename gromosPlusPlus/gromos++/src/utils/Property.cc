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

// 	$Id$	

//---Property Class-----------------------------------
#include "Property.h"

#include <cassert>
#include <map>
#include <math.h>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <cstdio>

#include "Value.h"
#include "AtomSpecifier.h"
#include "CremerPople.h"
#include "ExpressionParser.h"
#include "../gmath/Vec.h"
#include "../gcore/Bond.h"
#include "../gcore/Angle.h"
#include "../gcore/Dihedral.h"
#include "../utils/AtomSpecifier.h"
#include "../bound/Boundary.h"
#include "../gmath/Stat.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../fit/PositionUtils.h"

using namespace gcore;
using namespace std;
using namespace utils;

namespace utils {

  Property::Property(gcore::System &sys, bound::Boundary * pbc)
  : d_type("Property"),
  d_atom(sys),
  d_sys(&sys),
  d_pbc(pbc) {
  }

  Property::~Property() {
  }

  void Property::parse(std::vector<std::string> const & arguments, int x) {
    if (int(arguments.size()) < REQUIREDARGUMENTS)
      throw Exception(" too few arguments.\n");

    if (REQUIREDARGUMENTS > 0)
      parseAtoms(arguments[0], x);

    d_arg.resize(arguments.size() - 1, Value(0.0));

    for (unsigned int i = 1; i < arguments.size(); ++i)
      d_arg[i - 1].parse(arguments[i]);
  }

  void Property::parse(AtomSpecifier const & atmspc) {
    d_atom = atmspc;
  }

  void Property::parseAtoms(std::string atoms, int x) {
    d_atom.addSpecifier(atoms, x);
  }

  std::string Property::toTitle()const {
    std::ostringstream os;
    os << d_type << "%" << d_atom.toString()[0];
    return os.str();
  }

  std::ostream & operator<<(std::ostream &os, Property const & p) {
    os << p.toString();
    return os;
  }

  int Property::getTopologyType(gcore::System const &sys) {
    // standard implementation
    // check whether all atoms belong to the same molecule
    if (atoms().size()) {
      int the_mol = atoms().mol(0);

      for (unsigned int m = 1; m < atoms().size(); ++m)
        if (atoms().mol(m) != the_mol) return -1;

      return findTopologyType(sys.mol(the_mol).topology());
    } else return -1;
  }

  int Property::findTopologyType(gcore::MoleculeTopology const &mol_topo) {
    return -1;
  }

  std::string Property::average()const {
    std::ostringstream os;
    os.precision(8);
    os.setf(std::ios::fixed, std::ios::floatfield);

    if (d_scalar_stat.n()) {
      os << std::setw(15) << d_scalar_stat.ave()
              << std::setw(15) << d_scalar_stat.rmsd()
              << std::setw(15) << d_scalar_stat.ee();
    }
    if (d_vector_stat.n()) {
      // std::cerr << "getting vector stat: ave" << std::endl;
      os << std::setw(15) << gmath::v2s(d_vector_stat.ave());
      // std::cerr << "getting vector stat: rmsd" << std::endl;      
      os << std::setw(15) << gmath::v2s(d_vector_stat.rmsd());
      // std::cerr << "getting vector stat: ee" << std::endl;      
      // os << std::setw(15) << gmath::v2s(d_vector_stat.ee());
    }

    return os.str();
  }

  void Property::addValue(Value const & v) {
    switch (v.type()) {
      case val_scalar:
        d_scalar_stat.addval(v.scalar());
        break;
      case val_vector:
      case val_vecspec:
        d_vector_stat.addval(v.vec());
        break;
      default:
        throw Exception("bad value type");
    }
  }

  Value Property::nearestImageDistance(const Value & first, const Value & second) const{
    return second - first;
  }

  //---AverageProperty Class------------------------------------

  AverageProperty::AverageProperty
  (
          gcore::System &sys,
          bound::Boundary * pbc)
  : Property(sys, pbc) {
    d_type = "Average";
    REQUIREDARGUMENTS = 0;
  }

  Value const & AverageProperty::calc() {
    // empty
    d_single_scalar_stat = gmath::Stat<double>();
    d_single_vector_stat = gmath::Stat<gmath::Vec > ();

    for (unsigned int i = 0; i < d_property.size(); ++i) {
      Value const & v = d_property[i]->calc();

      switch (v.type()) {
        case val_scalar:
          d_single_scalar_stat.addval(v.scalar());
          break;
        case val_vector:
        case val_vecspec:
          d_single_vector_stat.addval(v.vec());
          break;
        default:
          throw Exception("wrong value type");
      }
    }

    if (d_single_scalar_stat.n())
      d_scalar_stat.addval(d_single_scalar_stat.ave());
    if (d_single_vector_stat.n())
      d_vector_stat.addval(d_single_vector_stat.ave());

    if (d_single_scalar_stat.n() < d_single_vector_stat.n())
      d_value = d_single_scalar_stat.ave();
    else
      d_value = d_single_vector_stat.ave();

    return d_value;
  }

  std::string AverageProperty::toTitle()const {
    std::ostringstream os;
    os << "<";
    bool first = true;
    for (unsigned int i = 0; i < d_property.size(); ++i) {
      if (!first)
        os << ",";
      first = false;
      os << d_property[i]->toTitle();
    }
    os << ">";
    return os.str();
  }

  std::string AverageProperty::toString()const {
    std::ostringstream os;

    if (d_single_scalar_stat.n()) {

      Value av = Value(d_single_scalar_stat.ave());
      Value rmsd = Value(d_single_scalar_stat.rmsd());
      Value ee = Value(d_single_scalar_stat.ee());

      os << av.toString() << " " << rmsd.toString() << " " << ee.toString();

      if (d_single_vector_stat.n()) os << "\t";
    }

    if (d_single_vector_stat.n()) {

      Value av = Value(d_single_vector_stat.ave());
      Value rmsd = Value(d_single_vector_stat.rmsd());
      Value ee = Value(d_single_vector_stat.ee());

      os << av.toString() << " " << rmsd.toString() << " " << ee.toString();
    }
    return os.str();
  }

  //---DistanceProperty Class------------------------------------

  DistanceProperty::DistanceProperty
  (
          gcore::System &sys,
          bound::Boundary * pbc)
  : Property(sys, pbc) {
    d_type = "Distance";
    REQUIREDARGUMENTS = 1;
  }

  void DistanceProperty::parse(std::vector<std::string> const & arguments, int x) {
    Property::parse(arguments, x);

    // for a distance, we should just have two atoms here...
    if (atoms().size() != 2) {
      std::cerr << "arguments:\n";
      for (unsigned int i = 0; i < arguments.size(); ++i)
        std::cerr << "\t" << arguments[i] << std::endl;
      std::cerr << "x = " << x << "\n" << std::endl;
      std::cerr << "resulted in a atom size of " << atoms().size() << std::endl;
      std::cerr << atoms().toString()[0] << std::endl;

      throw Exception("wrong number of atoms for a distance.\n");
    }
  }

  void DistanceProperty::parse(AtomSpecifier const & atmspc) {
    // for a distance, we should just have to atoms here...
    if (atmspc.size() != 2)
      throw Exception("wrong number of atoms for a distance.\n");
    Property::parse(atmspc);
  }

  Value const & DistanceProperty::calc() {
    gmath::Vec tmp = atoms().pos(0) -
            d_pbc->nearestImage(atoms().pos(0), atoms().pos(1), d_sys->box());

    const double d = tmp.abs();

    d_value = d;
    addValue(d_value);

    return d_value;
  }

  int DistanceProperty::findTopologyType
  (
          gcore::MoleculeTopology const &mol_topo
          ) {
    int a, b;

    if (atoms().atom(0) < atoms().atom(1)) {
      a = atoms().atom(0);
      b = atoms().atom(1);
    } else {
      a = atoms().atom(1);
      b = atoms().atom(0);
    }
    BondIterator bi(mol_topo);
    while (bi)
      if (bi()[0] == a && bi()[1] == b) return bi().type();
      else ++bi;

    return -1;
  }

  //---AngleProperty Class------------------------------------------------------------

  AngleProperty::AngleProperty(gcore::System &sys, bound::Boundary * pbc)
  : Property(sys, pbc) {
    d_type = "Angle";
    REQUIREDARGUMENTS = 1;
  }

  void AngleProperty::parse(std::vector<std::string> const & arguments, int x) {
    Property::parse(arguments, x);

    // it's an angle, therefore 3 atoms
    if (atoms().size() != 3)
      throw Exception("wrong number of atoms for an angle.\n");
  }

  void AngleProperty::parse(AtomSpecifier const & atmspc) {
    // it's an angle, therefore 3 atoms
    if (atmspc.size() != 3)
      throw Exception("wrong number of atoms for an angle.\n");
    Property::parse(atmspc);
  }

  Value const & AngleProperty::calc() {
    gmath::Vec tmpA = atoms().pos(0)
            - d_pbc->nearestImage(atoms().pos(0), atoms().pos(1), d_sys->box());
    gmath::Vec tmpB = atoms().pos(2)
            - d_pbc->nearestImage(atoms().pos(2), atoms().pos(1), d_sys->box());

    const double d = acos((tmpA.dot(tmpB)) / (tmpA.abs() * tmpB.abs()))*180 / M_PI;

    d_value = d;
    addValue(d_value);

    return d_value;
  }

  int AngleProperty::findTopologyType(gcore::MoleculeTopology const &mol_topo) {
    int a, b, c;

    if (atoms().atom(0) < atoms().atom(2)) {
      a = atoms().atom(0);
      b = atoms().atom(1);
      c = atoms().atom(2);
    } else {
      a = atoms().atom(2);
      b = atoms().atom(1);
      c = atoms().atom(0);
    }
    AngleIterator ai(mol_topo);
    while (ai)
      if (ai()[0] == a && ai()[1] == b && ai()[2] == c) return ai().type();
      else ++ai;

    return -1;
  }

  Value AngleProperty::nearestImageDistance(const Value & first, const Value & second) const {
    double diff = second.scalar() - first.scalar();
    while (diff >= 180.0) diff -= 360.0;
    while (diff < -180.0) diff += 360.0;
    return Value(diff);
  }

  //---TorsionProperty Class------------------------------------

  TorsionProperty::TorsionProperty(gcore::System &sys, bound::Boundary * pbc) :
  Property(sys, pbc) {
    d_type = "Torsion";
    REQUIREDARGUMENTS = 1;
  }

  TorsionProperty::~TorsionProperty() {
  }

  void TorsionProperty::parse(std::vector<std::string> const & arguments, int x) {
    Property::parse(arguments, x);

    // it's a torsion, therefore 4 atoms needed
    if (d_atom.size() != 4)
      throw Exception("wrong number of atoms for torsion.\n");
  }

  void TorsionProperty::parse(AtomSpecifier const & atmspc) {
    // it's a torsion, therefore 4 atoms needed
    if (atmspc.size() != 4)
      throw Exception("wrong number of atoms for torsion.\n");
    Property::parse(atmspc);
  }

  Value const & TorsionProperty::calc() {
    gmath::Vec tmpA = atoms().pos(0) - d_pbc->nearestImage(atoms().pos(0), atoms().pos(1), d_sys->box());
    gmath::Vec tmpB = atoms().pos(3) - d_pbc->nearestImage(atoms().pos(3), atoms().pos(2), d_sys->box());
    gmath::Vec tmpC = atoms().pos(2) - d_pbc->nearestImage(atoms().pos(2), atoms().pos(1), d_sys->box());

    gmath::Vec p1 = tmpA.cross(tmpC);
    gmath::Vec p2 = tmpB.cross(tmpC);

    double cosphi = ((p1.dot(p2)) / (p1.abs() * p2.abs()));

    if (cosphi > 1.0) cosphi = 1.0;
    if (cosphi <-1.0) cosphi = -1.0;
    
    double d = acos(cosphi)*180 / M_PI;
    
    gmath::Vec p3 = p1.cross(p2);
    if (p3.dot(tmpC) < 0) {
      d *= (-1);
    }
    
    if(d_scalar_stat.n() > 0) { // not the first caluclation => transformation needed
      double lastValue = d_scalar_stat.val(d_scalar_stat.n() - 1);
      double x = (d - lastValue) / 360.0;
      // get the nearest integer of x
      int ix;
      if(x > 0) {
        ix = int(x + 0.5);
      } else {
        ix = int(x - 0.5);
      }
      d = d - ix * 360.0;
    }
    
    d_value = d;
    addValue(d_value);
    // end of new angle definition

    return d_value;
  }

  int TorsionProperty::findTopologyType(gcore::MoleculeTopology const &mol_topo) {
    int a, b, c, d;

    if (atoms().atom(1) < atoms().atom(2)) {
      a = atoms().atom(0);
      b = atoms().atom(1);
      c = atoms().atom(2);
      d = atoms().atom(3);
    } else {
      a = atoms().atom(3);
      b = atoms().atom(2);
      c = atoms().atom(1);
      d = atoms().atom(0);
    }

    // assuming it is a dihedral...
    DihedralIterator di(mol_topo);
    while (di)
      if (di()[0] == a && di()[1] == b && di()[2] == c && di()[3] == d)
        return di().type();
      else ++di;

    return -1;
  }

  Value TorsionProperty::nearestImageDistance(const Value & first, const Value & second) const {
    double diff = second.scalar() - first.scalar();
    while (diff >= 180.0) diff -= 360.0;
    while (diff < -180.0) diff += 360.0;
    return Value(diff);
  }

  //---PeriodicTorsionProperty Class------------------------------------

  PeriodicTorsionProperty::PeriodicTorsionProperty(gcore::System &sys, bound::Boundary * pbc) :
  Property(sys, pbc) {
    d_type = "PeriodicTorsion";
    REQUIREDARGUMENTS = 1;
  }

  PeriodicTorsionProperty::~PeriodicTorsionProperty() {
  }

  void PeriodicTorsionProperty::parse(std::vector<std::string> const & arguments, int x) {
    Property::parse(arguments, x);

    // it's a torsion, therefore 4 atoms needed
    if (d_atom.size() != 4)
      throw Exception("wrong number of atoms for torsion.\n");
  }

  void PeriodicTorsionProperty::parse(AtomSpecifier const & atmspc) {
    // it's a torsion, therefore 4 atoms needed
    if (atmspc.size() != 4)
      throw Exception("wrong number of atoms for torsion.\n");
    Property::parse(atmspc);
  }

  Value const & PeriodicTorsionProperty::calc() {
    gmath::Vec tmpA = atoms().pos(0) - d_pbc->nearestImage(atoms().pos(0), atoms().pos(1), d_sys->box());
    gmath::Vec tmpB = atoms().pos(3) - d_pbc->nearestImage(atoms().pos(3), atoms().pos(2), d_sys->box());
    gmath::Vec tmpC = atoms().pos(2) - d_pbc->nearestImage(atoms().pos(2), atoms().pos(1), d_sys->box());

    gmath::Vec p1 = tmpA.cross(tmpC);
    gmath::Vec p2 = tmpB.cross(tmpC);

    double cosphi = ((p1.dot(p2)) / (p1.abs() * p2.abs()));

    if (cosphi > 1.0) cosphi = 1.0;
    if (cosphi <-1.0) cosphi = -1.0;
    
    double d = acos(cosphi)*180 / M_PI;
    
    // this is to know if the angle is positive or negative: arccos(cos(x)): [-1, 1] -> [0, pi], so the
    // sign has to be restored somehow
    gmath::Vec p3 = p1.cross(p2);
    if (p3.dot(tmpC) < 0) {
      d *= (-1);
    }
    
    if (d_arg.size()) {
      // shift into periodic range around zero value
      while (d < d_arg[0].scalar() - 180) {
        d += 360;
      }
      while (d > d_arg[0].scalar() + 180) {
        d -= 360;
      }
    }
    
    d_value = d;
    addValue(d_value);
    // end of new angle definition

    return d_value;
  }

  int PeriodicTorsionProperty::findTopologyType(gcore::MoleculeTopology const &mol_topo) {
    int a, b, c, d;

    if (atoms().atom(1) < atoms().atom(2)) {
      a = atoms().atom(0);
      b = atoms().atom(1);
      c = atoms().atom(2);
      d = atoms().atom(3);
    } else {
      a = atoms().atom(3);
      b = atoms().atom(2);
      c = atoms().atom(1);
      d = atoms().atom(0);
    }

    // assuming it is a dihedral...
    DihedralIterator di(mol_topo);
    while (di)
      if (di()[0] == a && di()[1] == b && di()[2] == c && di()[3] == d)
        return di().type();
      else ++di;

    return -1;
  }

  Value PeriodicTorsionProperty::nearestImageDistance(const Value & first, const Value & second) const {
    double diff = second.scalar() - first.scalar();
    while (diff >= 180.0) diff -= 360.0;
    while (diff < -180.0) diff += 360.0;
    return Value(diff);
  }
  
  //---CrossTorsionProperty Class------------------------------------

  CrossTorsionProperty::CrossTorsionProperty(gcore::System &sys, bound::Boundary * pbc) :
  Property(sys, pbc), d_atom2(sys) {
    d_type = "CrossTorsion";
    REQUIREDARGUMENTS = 2;
  }

  CrossTorsionProperty::~CrossTorsionProperty() {
  }

  void CrossTorsionProperty::parse(std::vector<std::string> const & arguments, int x) {
    if (int(arguments.size()) < REQUIREDARGUMENTS)
      throw Exception(" too few arguments.\n");

    parseAtoms(arguments[0], x);

    d_atom2.addSpecifier(arguments[1], x);

    d_arg.resize(arguments.size() - 2, Value(0.0));

    for (unsigned int i = 2; i < arguments.size(); ++i)
      d_arg[i - 1].parse(arguments[i]);

    // it's a cross torsion, therefore 8 atoms needed
    if (d_atom.size() != 4 || d_atom2.size() != 4)
      throw Exception("wrong number of atoms for cross torsion.\n");
  }

  void CrossTorsionProperty::parse(AtomSpecifier const & atmspc) {
    throw Exception("not implemented for cross torsion!\n");
  }

  Value const & CrossTorsionProperty::calc() {
    gmath::Vec tmpA = atoms().pos(0) - d_pbc->nearestImage(atoms().pos(0), atoms().pos(1), d_sys->box());
    gmath::Vec tmpB = atoms().pos(3) - d_pbc->nearestImage(atoms().pos(3), atoms().pos(2), d_sys->box());
    gmath::Vec tmpC = atoms().pos(2) - d_pbc->nearestImage(atoms().pos(2), atoms().pos(1), d_sys->box());

    gmath::Vec p1 = tmpA.cross(tmpC);
    gmath::Vec p2 = tmpB.cross(tmpC);

    double cosphi = ((p1.dot(p2)) / (p1.abs() * p2.abs()));

    if (cosphi > 1.0) cosphi = 1.0;
    if (cosphi <-1.0) cosphi = -1.0;
    double d1 = acos(cosphi)*180 / M_PI;

    gmath::Vec p3 = p1.cross(p2);
    if (p3.dot(tmpC) < 0)
      d1 = 360 - d1;

    if (d_arg.size()) {
      // shift into periodic range around zero value
      while (d1 < d_arg[0].scalar() - 180) {
        d1 += 360;
      }
      while (d1 > d_arg[0].scalar() + 180) {
        d1 -= 360;
      }
    }

    tmpA = atoms2().pos(0) - d_pbc->nearestImage(atoms2().pos(0), atoms2().pos(1), d_sys->box());
    tmpB = atoms2().pos(3) - d_pbc->nearestImage(atoms2().pos(3), atoms2().pos(2), d_sys->box());
    tmpC = atoms2().pos(2) - d_pbc->nearestImage(atoms2().pos(2), atoms2().pos(1), d_sys->box());

    p1 = tmpA.cross(tmpC);
    p2 = tmpB.cross(tmpC);

    cosphi = ((p1.dot(p2)) / (p1.abs() * p2.abs()));

    if (cosphi > 1.0) cosphi = 1.0;
    if (cosphi <-1.0) cosphi = -1.0;
    double d2 = acos(cosphi)*180 / M_PI;

    p3 = p1.cross(p2);
    if (p3.dot(tmpC) < 0)
      d2 = 360 - d2;

    if (d_arg.size()) {
      // shift into periodic range around zero value
      while (d2 < d_arg[0].scalar() - 180) {
        d2 += 360;
      }
      while (d2 > d_arg[0].scalar() + 180) {
        d2 -= 360;
      }
    }

    d_value = d1 + d2;
    addValue(d_value);

    return d_value;
  }

  int CrossTorsionProperty::findTopologyType(gcore::MoleculeTopology const &mol_topo) {
    // not implemented!
    return -1;
  }

  std::string CrossTorsionProperty::toTitle()const {
    std::ostringstream os;
    os << d_type << "%" << d_atom.toString()[0] << "%" << d_atom2.toString()[0];
    return os.str();
  }

  Value CrossTorsionProperty::nearestImageDistance(const Value & first, const Value & second) const {
    double diff = second.scalar() - first.scalar();
    while (diff >= 180.0) diff -= 360.0;
    while (diff < -180.0) diff += 360.0;
    return Value(diff);
  }

  //---JValueProperty Class------------------------------------

  JValueProperty::JValueProperty(gcore::System &sys, bound::Boundary * pbc) :
  Property(sys, pbc) {
    d_type = "JValue";
    REQUIREDARGUMENTS = 1;
  }

  JValueProperty::~JValueProperty() {
  }

  void JValueProperty::parse(std::vector<std::string> const & arguments, int x) {
    Property::parse(arguments, x);

    // it's a torsion, therefore 4 atoms needed
    if (d_atom.size() != 4)
      throw Exception("wrong number of atoms for j-value: 4 needed.\n");

    if (d_arg.size() == 0) {
      d_arg.resize(4, Value(0.0));
      // take default parameters
      d_arg[0].parse("6.4"); // a
      d_arg[1].parse("-1.4"); // b
      d_arg[2].parse("1.9"); // c
      d_arg[3].parse("0.0"); // delta
    } else if (d_arg.size() != 4) {
      throw Exception("not four arguments for J value.\n");
    }
  }

  Value const & JValueProperty::calc() {
    gmath::Vec tmpA = atoms().pos(0) - d_pbc->nearestImage(atoms().pos(0), atoms().pos(1), d_sys->box());
    gmath::Vec tmpB = atoms().pos(3) - d_pbc->nearestImage(atoms().pos(3), atoms().pos(2), d_sys->box());
    gmath::Vec tmpC = atoms().pos(2) - d_pbc->nearestImage(atoms().pos(2), atoms().pos(1), d_sys->box());

    gmath::Vec p1 = tmpA.cross(tmpC);
    gmath::Vec p2 = tmpB.cross(tmpC);

    double cosphi = ((p1.dot(p2)) / (p1.abs() * p2.abs()));

    if (cosphi > 1.0) cosphi = 1.0;
    if (cosphi <-1.0) cosphi = -1.0;
    double d = acos(cosphi)*180 / M_PI;

    gmath::Vec p3 = p1.cross(p2);
    if (p3.dot(tmpC) < 0)
      d = 360 - d;

    double a = d_arg[0].scalar();
    double b = d_arg[1].scalar();
    double c = d_arg[2].scalar();
    double delta = (d_arg[3].scalar() / 180) * M_PI;
    //std::cerr << "a: " << a << "\t b: " << b << "\tc: " <<c << std::endl;

    double phi = (d / 180) * M_PI;
    //std::cerr << "d: " << d << "\t phi: " << phi << std::endl;
    const double the_cos = cos(phi + delta);
    double J = a * the_cos * the_cos + b * the_cos + c;
    //std::cerr << "J: " << J << std::endl;

    d_value = J;
    addValue(d_value);

    return d_value;
  }

  //---VectorOrderProperty Class------------------------------------------------------------

  VectorOrderProperty::VectorOrderProperty(gcore::System &sys, bound::Boundary * pbc) :
  Property(sys, pbc),
  d_vec1(sys, pbc),
  d_vec2(sys, pbc) {
    // the first one we read ourselves...
    REQUIREDARGUMENTS = 2;
    d_type = "VectorOrder";
  }

  VectorOrderProperty::~VectorOrderProperty() {
  }

  void VectorOrderProperty::parse(std::vector<std::string> const & arguments, int x) {
    if (arguments.size() != 2) {
      throw Exception(" VectorOrder needs two vector specifiers as argument.\n");
    }

    std::map<std::string, utils::Value> vars;
    if (x != -1) {
      vars["x"] = Value(x);
    }

    d_vec1.setSpecifier(arguments[0], vars);
    d_vec2.setSpecifier(arguments[1], vars);
  }

  Value const & VectorOrderProperty::calc() {
    gmath::Vec tmpA = d_vec2();
    gmath::Vec d_axis = d_vec1();

    const double d = acos((tmpA.dot(d_axis)) / (tmpA.abs() * d_axis.abs()))*180 / M_PI;

    d_value = d;
    addValue(d_value);

    return d_value;
  }

  std::string VectorOrderProperty::toTitle()const {

    std::ostringstream s;
    s << "o%"
            << d_vec1.toString() << "%" << d_vec2.toString();

    return s.str();
  }

  int VectorOrderProperty::findTopologyType(gcore::MoleculeTopology const &mol_topo) {
    // not implemented
    return -1;
  }

  //---VectorOrderParamProperty Class------------------------------------------------------------

  VectorOrderParamProperty::VectorOrderParamProperty(gcore::System &sys, bound::Boundary * pbc) :
  Property(sys, pbc),
  d_vec1(sys, pbc),
  d_vec2(sys, pbc) {
    // the first one we read ourselves...
    REQUIREDARGUMENTS = 2;
    d_type = "VectorOrderParam";
  }

  VectorOrderParamProperty::~VectorOrderParamProperty() {
  }

  void VectorOrderParamProperty::parse(std::vector<std::string> const & arguments, int x) {
    if (arguments.size() != 2) {
      throw Exception(" VectorOrderParam needs two vector specifiers as argument.\n");
    }

    std::map<std::string, utils::Value> vars;
    if (x != -1) {
      vars["x"] = Value(x);
    }

    d_vec1.setSpecifier(arguments[0], vars);
    d_vec2.setSpecifier(arguments[1], vars);
  }

  Value const & VectorOrderParamProperty::calc() {
    gmath::Vec tmpA = d_vec2();
    gmath::Vec d_axis = d_vec1();

    const double cosa = tmpA.dot(d_axis) / (tmpA.abs() * d_axis.abs());
    const double d = 0.5 * (3 * cosa * cosa - 1);

    d_value = d;
    addValue(d_value);

    return d_value;
  }

  std::string VectorOrderParamProperty::toTitle()const {

    std::ostringstream s;
    s << "op%"
            << d_vec1.toString() << "%" << d_vec2.toString();

    return s.str();
  }

  int VectorOrderParamProperty::findTopologyType(gcore::MoleculeTopology const &mol_topo) {
    // not implemented
    return -1;
  }

  //---PseudoRotation Class------------------------------------

  PseudoRotationProperty::PseudoRotationProperty(gcore::System &sys, bound::Boundary * pbc) :
  Property(sys, pbc),
  d_sin36sin72(::sin(M_PI / 5.0) + ::sin(2.0 * M_PI / 5.0)) {
    d_type = "PseudoRotation";
    REQUIREDARGUMENTS = 1;
  }

  PseudoRotationProperty::~PseudoRotationProperty() {
  }

  void PseudoRotationProperty::parse(std::vector<std::string> const & arguments, int x) {
    Property::parse(arguments, x);

    // it's a pseudo rotation, therefore 5 atoms needed
    if (d_atom.size() != 5)
      throw Exception(" wrong number of atoms for pseudo rotation.\n");
  }

  Value const & PseudoRotationProperty::calc() {
    // first calculate the four dihedrals
    double torsion[5];
    torsion[0] = _calcDihedral(0, 1, 2, 3);
    torsion[1] = _calcDihedral(1, 2, 3, 4);
    torsion[2] = _calcDihedral(2, 3, 4, 0);
    torsion[3] = _calcDihedral(3, 4, 0, 1);
    torsion[4] = _calcDihedral(4, 0, 1, 2);

    for (int i = 0; i < 5; ++i)
      if (torsion[i] > 180.0) torsion[i] -= 360.0;

    double factor = (torsion[2] + torsion[4]) - (torsion[1] + torsion[3]);

    factor /= (2.0 * torsion[0] * d_sin36sin72);

    double d = atan(factor)*180 / M_PI;

    if (torsion[0] < 0.0) d += 180;

    d_value = d;
    addValue(d_value);

    return d_value;
  }

  double PseudoRotationProperty::
  _calcDihedral(int const a, int const b, int const c, int const d) {
    gmath::Vec tmpA = atoms().pos(a) - d_pbc->nearestImage(atoms().pos(a), atoms().pos(b), d_sys->box());
    gmath::Vec tmpB = atoms().pos(d) - d_pbc->nearestImage(atoms().pos(d), atoms().pos(c), d_sys->box());
    gmath::Vec tmpC = atoms().pos(c) - d_pbc->nearestImage(atoms().pos(c), atoms().pos(b), d_sys->box());

    gmath::Vec p1 = tmpA.cross(tmpC);
    gmath::Vec p2 = tmpB.cross(tmpC);

    double cosphi = ((p1.dot(p2)) / (p1.abs() * p2.abs()));

    double value = acos(cosphi)*180 / M_PI;

    gmath::Vec p3 = p1.cross(p2);
    if (p3.dot(tmpC) < 0)
      value = 360 - value;

    return value;
  }

  int PseudoRotationProperty::findTopologyType(gcore::MoleculeTopology const &mol_topo) {
    // maybe give the residue number?
    if (mol_topo.resNum(atoms().atom(0)) == mol_topo.resNum(atoms().atom(1)) &&
            mol_topo.resNum(atoms().atom(0)) == mol_topo.resNum(atoms().atom(2)) &&
            mol_topo.resNum(atoms().atom(0)) == mol_topo.resNum(atoms().atom(3)) &&
            mol_topo.resNum(atoms().atom(0)) == mol_topo.resNum(atoms().atom(4)))
      return mol_topo.resNum(atoms().atom(0));
    else
      return -1;
  }

  //---PuckerAmplitude Class------------------------------------

  PuckerAmplitudeProperty::PuckerAmplitudeProperty(gcore::System &sys, bound::Boundary * pbc) :
  PseudoRotationProperty(sys, pbc) {
    d_type = "PuckerAmplitude";
    //REQUIREDARGUMENTS = 1;
    //d_sin36sin72 = sin(M_PI/5.0) + sin(2.0*M_PI/5.0);
  }

  PuckerAmplitudeProperty::~PuckerAmplitudeProperty() {
  }

  Value const & PuckerAmplitudeProperty::calc() {
    // first calculate the four dihedrals
    double torsion[5];
    torsion[0] = _calcDihedral(0, 1, 2, 3);
    torsion[1] = _calcDihedral(1, 2, 3, 4);
    torsion[2] = _calcDihedral(2, 3, 4, 0);
    torsion[3] = _calcDihedral(3, 4, 0, 1);
    torsion[4] = _calcDihedral(4, 0, 1, 2);
    for (int i = 0; i < 5; ++i)
      if (torsion[i] > 180.0) torsion[i] -= 360.0;

    double factor = (torsion[2] + torsion[4]) - (torsion[1] + torsion[3]);
    factor /= (2.0 * torsion[0] * d_sin36sin72);

    double pr = atan(factor);
    double t0 = torsion[0];

    if (torsion[0] < 0) {
      pr += M_PI;
    }

    const double d = t0 / ::cos(pr);

    d_value = d;
    addValue(d_value);

    return d_value;
  }

  //---CPAmplitude Coordinates Class------------------------------------

  CPAmplitudeProperty::CPAmplitudeProperty(gcore::System &sys, bound::Boundary * pbc) :
  Property(sys, pbc) {
    d_type = "CPAmplitude";
    REQUIREDARGUMENTS = 1;
  }

  CPAmplitudeProperty::~CPAmplitudeProperty() {
  }
  
  void CPAmplitudeProperty::parse(std::vector<std::string> const & arguments, int x) {
    Property::parse(arguments, x);

    //cremer pople coordinates with 6 atoms
    if (d_atom.size() != 6)
      throw Exception(" wrong number of atoms for pseudo rotation.\n");
  }
  
  Value const & CPAmplitudeProperty::calc() {  
    std::vector<gmath::Vec> R(6); //atom coordinates in new coordinate system, 
    gmath::Vec Rp, Rdp, RpxRdp; //R', R''
    gmath::Vec Rcm, temp; //center of geometry, normal vector, temporary vector
  
    R[0]=atoms().pos(0);

    for (unsigned int i = 1; i < 6; ++i) {
      R[i]=d_pbc->nearestImage(R[i-1], atoms().pos(i), d_sys->box());
    } 
    util::cremerpople::calcCOM(R, Rcm);
    util::cremerpople::shiftToCOM(R, Rcm);
    util::cremerpople::calcPlane(R, Rp, Rdp, RpxRdp);
    
    //displacement vector in z direction and total puckering amplitude Q
    std::vector<double> zeta(6);
    d_value = util::cremerpople::calcZeta(R, RpxRdp, zeta);
    addValue(d_value);

    return d_value;
  }

  //---CPTheta Coordinates Class------------------------------------

  CPThetaProperty::CPThetaProperty(gcore::System &sys, bound::Boundary * pbc) :
  Property(sys, pbc) {
    d_type = "CPTheta";
    REQUIREDARGUMENTS = 1;
  }

  CPThetaProperty::~CPThetaProperty() {
  }
  
  void CPThetaProperty::parse(std::vector<std::string> const & arguments, int x) {
    Property::parse(arguments, x);

    //cremer pople coordinates with 6 atoms
    if (d_atom.size() != 6)
      throw Exception(" wrong number of atoms for pseudo rotation.\n");
  }
  
  Value const & CPThetaProperty::calc() {  
    std::vector<gmath::Vec> R(6); //atom coordinates in new coordinate system, 
    gmath::Vec Rp, Rdp, RpxRdp; //R', R''
    gmath::Vec Rcm, temp; //center of geometry, normal vector, temporary vector
  
    R[0]=atoms().pos(0);

    for (unsigned int i = 1; i < 6; ++i) {
      R[i]=d_pbc->nearestImage(R[i-1], atoms().pos(i), d_sys->box());
    } 
    util::cremerpople::calcCOM(R, Rcm);
    util::cremerpople::shiftToCOM(R, Rcm);
    util::cremerpople::calcPlane(R, Rp, Rdp, RpxRdp);
    
    //displacement vector in z direction and total puckering amplitude Q
    std::vector<double> zeta(6);
    util::cremerpople::calcZeta(R, RpxRdp, zeta);
  
    //theta
    double theta = util::cremerpople::calcTheta(zeta);
    d_value = theta/M_PI*180.0;
    addValue(d_value);

    return d_value;
  }

  //---CPPhi Coordinates Class------------------------------------

  CPPhiProperty::CPPhiProperty(gcore::System &sys, bound::Boundary * pbc) :
  Property(sys, pbc) {
    d_type = "CPPhi";
    REQUIREDARGUMENTS = 1;
  }

  CPPhiProperty::~CPPhiProperty() {
  }
  
  void CPPhiProperty::parse(std::vector<std::string> const & arguments, int x) {
    Property::parse(arguments, x);

    //cremer pople coordinates with 6 atoms
    if (d_atom.size() != 6)
      throw Exception(" wrong number of atoms for pseudo rotation.\n");
  }
  
  Value const & CPPhiProperty::calc() {  
    std::vector<gmath::Vec> R(6); //atom coordinates in new coordinate system, 
    gmath::Vec Rp, Rdp, RpxRdp; //R', R''
    gmath::Vec Rcm, temp; //center of geometry, normal vector, temporary vector
  
    R[0]=atoms().pos(0);

    for (unsigned int i = 1; i < 6; ++i) {
      R[i]=d_pbc->nearestImage(R[i-1], atoms().pos(i), d_sys->box());
    } 
    util::cremerpople::calcCOM(R, Rcm);
    util::cremerpople::shiftToCOM(R, Rcm);
    util::cremerpople::calcPlane(R, Rp, Rdp, RpxRdp);

    //displacement vector in z direction and total puckering amplitude Q
    std::vector<double> zeta(6);
    util::cremerpople::calcZeta(R, RpxRdp, zeta);
  
    //phi
    double phi = util::cremerpople::calcPhi(zeta);
    d_value = phi/M_PI*180.0;
    addValue(d_value);

    return d_value;
  }

  //---ExpressionProperty Class------------------------------------------------------------

  ExpressionProperty::ExpressionProperty(gcore::System &sys, bound::Boundary * pbc) :
  Property(sys, pbc),
  d_parser(&sys, pbc) {
    // the first one we read ourselves...
    REQUIREDARGUMENTS = 1;
    d_type = "Expression";
  }

  ExpressionProperty::~ExpressionProperty() {
  }

  void ExpressionProperty::parse(std::vector<std::string> const & arguments, int x) {
    d_expr.clear();
    d_var["x"] = x;
    d_parser.parse_expression(arguments[0], d_var, d_expr);

    // nothing left to parse ...
    // Property::parse(count - 2, &arguments[1]);
  }

  Value const & ExpressionProperty::calc() {
    d_value = d_parser.calculate(d_expr, d_var);
    addValue(d_value);

    return d_value;
  }

  std::string ExpressionProperty::toTitle()const {
    std::ostringstream s;
    s << "expr";
    return s.str();
  }

  int ExpressionProperty::findTopologyType(gcore::MoleculeTopology const &mol_topo) {
    // not implemented
    return -1;
  }

  //---HBProperty Class------------------------------------------------------------------

  HBProperty::HBProperty(gcore::System &sys, bound::Boundary * pbc) :
  Property(sys, pbc),
  /**  HB distances (2 if 3 center)
   */
  d1_hb(sys, pbc), d2_hb(sys, pbc),
  /** HB angles (3 if 3 center)
   */
  a1_hb(sys, pbc), a2_hb(sys, pbc), a3_hb(sys, pbc),
  /** HB improper (1 if 3 center)
   */
  i1_hb(sys, pbc) {
    d_type = "HB";
    REQUIREDARGUMENTS = 1;
  }

  HBProperty::~HBProperty() {
  }

  void HBProperty::parse(std::vector<std::string> const & arguments, int x) {
    Property::parse(arguments, x);

    if (atoms().size() < 3 || atoms().size() > 4)
      throw Exception(" no enough of too many atoms for hydrogen bond determination.\n");

    // 'atoms' order should be D-H-A1(-A2)
    if (atoms().mass(1) != 1.00800)
      throw Exception(" the second atom must be a Hydrogen.\n");

    // we have a normal hydrogen bond
    if (atoms().size() == 3) {
      is3c = false;
      AtomSpecifier tmp(*atoms().sys()); // need to reorder

      // D-H-O angle
      tmp.addAtom(atoms().mol(0), atoms().atom(0));
      tmp.addAtom(atoms().mol(1), atoms().atom(1));
      tmp.addAtom(atoms().mol(2), atoms().atom(2));
      a1_hb.parse(tmp);

      // H-A distance
      tmp.clear();
      tmp.addAtom(atoms().mol(1), atoms().atom(1));
      tmp.addAtom(atoms().mol(2), atoms().atom(2));
      d1_hb.parse(tmp);

      if (d_arg.size() == 0) {
        d_arg.resize(2, Value(0.0));
        d_arg[0].parse("0.25"); //default values for hbond
        d_arg[1].parse("135"); // max dist && min ang
      } else if (d_arg.size() != 2) {
        throw Exception(" not two arguments for HB.\n");
      }
    } else if (atoms().size() == 4) { // 3 center hb
      is3c = true;
      AtomSpecifier tmp(*atoms().sys()); // need to reorder

      // torsion through planes
      tmp.addAtom(atoms().mol(1), atoms().atom(1));
      tmp.addAtom(atoms().mol(3), atoms().atom(3));
      tmp.addAtom(atoms().mol(2), atoms().atom(2));
      tmp.addAtom(atoms().mol(0), atoms().atom(0));
      i1_hb.parse(tmp);

      // H-A1 distance
      tmp.clear();
      tmp.addAtom(atoms().mol(1), atoms().atom(1));
      tmp.addAtom(atoms().mol(2), atoms().atom(2));
      d1_hb.parse(tmp);

      // H-A2 distance
      tmp.clear();
      tmp.addAtom(atoms().mol(1), atoms().atom(1));
      tmp.addAtom(atoms().mol(3), atoms().atom(3));
      d2_hb.parse(tmp);

      // D-H-A1 angle
      tmp.clear();
      tmp.addAtom(atoms().mol(0), atoms().atom(0));
      tmp.addAtom(atoms().mol(1), atoms().atom(1));
      tmp.addAtom(atoms().mol(2), atoms().atom(2));
      a1_hb.parse(tmp);

      // D-H-A2 angle
      tmp.clear();
      tmp.addAtom(atoms().mol(0), atoms().atom(0));
      tmp.addAtom(atoms().mol(1), atoms().atom(1));
      tmp.addAtom(atoms().mol(3), atoms().atom(3));
      a2_hb.parse(tmp);

      // A1-H-A2 angle
      tmp.clear();
      tmp.addAtom(atoms().mol(2), atoms().atom(2));
      tmp.addAtom(atoms().mol(1), atoms().atom(1));
      tmp.addAtom(atoms().mol(3), atoms().atom(3));
      a3_hb.parse(tmp);

      if (d_arg.size() == 0) {
        d_arg.resize(4, Value(0.0));
        d_arg[0].parse("0.27"); //default values for hbond
        d_arg[1].parse("90"); // max dist && min ang
        d_arg[2].parse("340"); // & min summ & max imp
        d_arg[3].parse("15");
      } else if (d_arg.size() != 4) {
        throw Exception(" not four arguments for three centered HB.\n");
      }
    }
  }

  Value const & HBProperty::calc() {
    if (is3c)
      if (d1_hb.calc().scalar() <= d_arg[0].scalar() &&
              d2_hb.calc().scalar() <= d_arg[0].scalar() &&
              a1_hb.calc().scalar() >= d_arg[1].scalar() &&
              a2_hb.calc().scalar() >= d_arg[1].scalar() &&
              (a1_hb.calc().scalar() +
              a2_hb.calc().scalar() +
              a3_hb.calc().scalar()) >= d_arg[2].scalar() &&
              i1_hb.calc().scalar() <= d_arg[3].scalar())
        d_value = 1;
      else
        d_value = 0;
    else // no 3c HB
      if (d1_hb.calc().scalar() <= d_arg[0].scalar() &&
            a1_hb.calc().scalar() >= d_arg[1].scalar())
      d_value = 1;
    else
      d_value = 0;
    addValue(d_value);

    return d_value;
  }

  int HBProperty::findTopologyType(gcore::MoleculeTopology const &mol_topo) {
    if (is3c)
      return 1;
    else
      return 0;
  }

  std::string HBProperty::toTitle()const {
    std::ostringstream os;
    if (is3c)
      os << d_type << "%" << atoms().toString()[0]
      << "%" << d_arg[0].scalar()
      << "%" << d_arg[1].scalar()
      << "%" << d_arg[2].scalar()
      << "%" << d_arg[3].scalar();
    else
      os << d_type << "%" << atoms().toString()[0]
      << "%" << d_arg[0].scalar()
      << "%" << d_arg[1].scalar();
    return os.str();
  }
  //---StackingProperty Class------------------------------------------------------------------

  StackingProperty::StackingProperty(gcore::System &sys, bound::Boundary * pbc) :
  Property(sys, pbc), d_atoms1(sys), d_atoms2(sys) {
    d_type = "st";
    REQUIREDARGUMENTS = 2;
  }

  StackingProperty::~StackingProperty() {
  }

  void StackingProperty::parse(std::vector<std::string> const & arguments, int x) {
    if (int(arguments.size()) != REQUIREDARGUMENTS &&
            int(arguments.size()) != REQUIREDARGUMENTS + 2) {
      throw Exception(" wrong number of arguments for stacking property.");
    }

    // add the first arguments to the first plane - with subst.
    atoms1().addSpecifier(arguments[0], x);
    if (atoms1().size() < 3) {
      ostringstream msg;
      msg << " less than 3 atoms for first plane:";
      for (int i = 0; atoms1().size(); ++i)
        msg << " " << atoms1().name(i);
      throw Exception(msg.str());
    }
    // add the plane to the atoms: this is not really needed but the user may
    // expect something useful if accessing atoms()
    for (unsigned int i = 0; i < atoms1().size(); ++i)
      atoms().addAtom(atoms1().mol(i), atoms1().atom(i));

    // add the second arguments to the second plane - with subst.
    atoms2().addSpecifier(arguments[1], x);
    if (atoms2().size() < 3) {
      ostringstream msg;
      msg << " less than 3 atoms for second plane:";
      for (int i = 0; atoms2().size(); ++i)
        msg << " " << atoms2().name(i);
      throw Exception(msg.str());
    }
    for (unsigned int i = 0; i < atoms2().size(); ++i)
      atoms().addAtom(atoms2().mol(i), atoms2().atom(i));

    if (int(arguments.size()) == REQUIREDARGUMENTS) {
      // we only got the planes: set default values
      d_arg.resize(2, Value(0.0));
      d_arg[0].parse("0.5");
      d_arg[1].parse("30");
    } else {
      d_arg.resize(2, Value(0.0));
      d_arg[0].parse(arguments[2]);
      d_arg[1].parse(arguments[3]);
    }
  }

  Value const & StackingProperty::calc() {
    d_value = 0;
    // using the first three atoms (a, b, c) in the specifier calculate the norm vector
    // to the plane through these three atoms:
    // i.e. a_b = a - b; a_c = a - c; but we gather to a to be on the safe side.
    const gmath::Vec a_b_1 = atoms1().pos(0) -
            d_pbc->nearestImage(atoms1().pos(0), atoms1().pos(1), d_sys->box());
    const gmath::Vec a_c_1 = atoms1().pos(0) -
            d_pbc->nearestImage(atoms1().pos(0), atoms1().pos(2), d_sys->box());
    const gmath::Vec norm1 = a_b_1.cross(a_c_1);

    const gmath::Vec a_b_2 = atoms2().pos(0) -
            d_pbc->nearestImage(atoms2().pos(0), atoms2().pos(1), d_sys->box());
    const gmath::Vec a_c_2 = atoms2().pos(0) -
            d_pbc->nearestImage(atoms2().pos(0), atoms2().pos(2), d_sys->box());
    const gmath::Vec norm2 = a_b_2.cross(a_c_2);

    // now we have the norms - let's calculate the angles.
    const double angle = fabs(acos(norm1.dot(norm2) / (norm1.abs() * norm2.abs()))) * 180.0 / M_PI;

    // check whether the angle is lower than the upper bound
    // abort if not - so we save the cog calculation.
    if (angle > d_arg[1].scalar()) {
      addValue(d_value);
      return d_value;
    }

    // calculate the centres of geometries of both planes.
    const gmath::Vec cog1 = fit::PositionUtils::cog(*d_sys, atoms1());
    const gmath::Vec cog2 = fit::PositionUtils::cog(*d_sys, atoms2());

    const double d = (cog1 - d_pbc->nearestImage(cog1, cog2, d_sys->box())).abs();
    // check the distance and abort
    if (d > d_arg[0].scalar()) {
      addValue(d_value);
      return d_value;
    }

    // they do stack
    d_value = 1;
    addValue(d_value);
    return d_value;
  }

  int StackingProperty::findTopologyType(gcore::MoleculeTopology const &mol_topo) {
    return 0;
  }

  std::string StackingProperty::toTitle()const {
    std::ostringstream os;
    os << d_type << "%" << atoms1().toString()[0]
            << "%" << atoms2().toString()[0]
            << "%" << d_arg[0].scalar()
            << "%" << d_arg[1].scalar();
    return os.str();
  }

} // utils


