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
#include "PropertyContainer.h"

#include <cassert>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "parse.h"
#include "../gcore/System.h"
#include "../bound/Boundary.h"

using namespace gcore;
using namespace std;
using namespace utils;

namespace utils
{
  //---PropertyContainer Class------------------------------------

  PropertyContainer::PropertyContainer() :
    vector<Property *>()
  {
    d_sys = NULL;
  }
  
  PropertyContainer::PropertyContainer(gcore::System &sys, bound::Boundary *pbc)
    : vector<Property *>(),
      d_sys(&sys),
      d_pbc(pbc)
  {
  }
  
   PropertyContainer::~PropertyContainer()
  {
    // delete all properties
    // (this is slightly undefined behaviour, as the deleted pointers
    //  are left in the vector (for clear))
    for(std::vector<Property *>::iterator it = this->begin(),
	  to = this->end(); it != to; ++it){
      delete *it;
    }
  }
  
  void PropertyContainer::reinitialize(gcore::System &sys, bound::Boundary *pbc)
  {
    for(std::vector<Property *>::iterator it = this->begin(),
	  to = this->end(); it != to; ++it){
      delete *it;
    }
    this->clear();
    d_sys=&sys;
    d_pbc=pbc;
  }
  
  int PropertyContainer::addSpecifier(std::string s)
  {
    // std::cerr << "adding: " << s << std::endl;
    parse(s);
    return size();
  }

  // split properties on ';'
  void PropertyContainer::parse(std::string s)
  {
    parse_multiple(s, *this);
  }

  void PropertyContainer::parse_multiple(std::string s, std::vector<Property *> & prop)
  {
    std::string::size_type it = find_par(s, ' ', 0);
    
    if (it == std::string::npos){
      if (s.find_first_not_of(" ") != std::string::npos)
	parse_single(s, prop);
    }
    else{
      if (s.find_first_not_of(" ") < it)
	      parse_single(s.substr(0, it), prop);
      parse_multiple(s.substr(it+1, std::string::npos), prop);
    }
  }

  // separate distribution metaproperty

  void PropertyContainer::parse_single(std::string s, std::vector<Property *> & prop) {
    // check if we have an average property
    if (s[0] == '<') {
      char bra = '<';
      std::string::size_type end_dist = find_matching_bracket(s, bra, 1);
      if (end_dist == std::string::npos) {
        throw Exception("AverageProperty: closing bracket wrong!");
      }

      parse_average(s.substr(1, end_dist - 2));
      // nothing can follow an average (as single property)
      return;
    }

    // look for a range definition for x
    std::vector<int> x_range;
    std::string::size_type var_start = s.find('|');
    if (var_start != std::string::npos) {

      std::string::size_type var_start2
          = s.substr(var_start + 1, std::string::npos).find('=');

      if (var_start2 == std::string::npos)
        throw Exception("could not parse variable definition!");

      std::string var_part = s.substr(var_start + var_start2 + 2,
          std::string::npos);
      s = s.substr(0, var_start);

      parse_range(var_part, x_range);
    }

    // extract type
    std::string type;
    std::string::size_type it = s.find('%');
    if (it == std::string::npos) {
      std::cerr << "property: " << s << std::endl;
      throw Exception
          (" invalid property-specifier.\nSyntax: <type>%<atomspecifier>[%...]\n");
    }

    type = s.substr(0, it);

    std::vector<std::string> arguments;

    // now separate the string on % into argument tokens
    for (int c = 0; c < Property::MAXARGUMENTS; ++c) {

      std::string::size_type last = it;

      it = s.find('%', last + 1);
      if (it == std::string::npos) {
        // the last argument
        arguments.push_back(s.substr(last + 1, std::string::npos));
        break;
      }

      arguments.push_back(s.substr(last + 1, it - last - 1));
    }

    if (x_range.size()) {
      for (unsigned int x = 0; x < x_range.size(); ++x)
        prop.push_back(createProperty(type, arguments, x_range[x]));
    } else {
      prop.push_back(createProperty(type, arguments, -1));
    }
  }
// end of parse_single

// create an average property

void PropertyContainer::parse_average(std::string s) {

  AverageProperty *av = new AverageProperty(*d_sys, d_pbc);
  push_back(av);

  parse_multiple(s, av->properties());

}

  ////////////////////////////////////////////////////////////
  // Create a Property
  ////////////////////////////////////////////////////////////
  Property * PropertyContainer::createProperty
  (
   std::string type,
   std::vector<std::string> const & arguments,
   int x
   )
  {
    // this method has to know all existing property types
    // when adding user properties, a child class of PropertyContainer
    // has to be added which extends this method for the new properties
    if (type == "d"){ 
      DistanceProperty *p = new DistanceProperty(*d_sys, d_pbc);
      p->parse(arguments, x);
      return p;
    }
    
    if (type == "a"){
      AngleProperty *p = new AngleProperty(*d_sys, d_pbc);
      p->parse(arguments, x);
      return p;
    }
    
    if (type == "t"){
      TorsionProperty *p = new TorsionProperty(*d_sys, d_pbc);
      p->parse(arguments, x);
      return p;
    }
    
    if (type == "tp") {
      PeriodicTorsionProperty *p = new PeriodicTorsionProperty(*d_sys, d_pbc);
      p->parse(arguments, x);
      return p;
    }
      
    if (type == "ct"){
      CrossTorsionProperty *p = new CrossTorsionProperty(*d_sys, d_pbc);
      p->parse(arguments, x);
      return p;
    }

    if (type == "j"){
      JValueProperty *p = new JValueProperty(*d_sys, d_pbc);
      p->parse(arguments, x);
      return p;
    }

    if (type == "o"){
      VectorOrderProperty *p = new VectorOrderProperty(*d_sys, d_pbc);
      p->parse(arguments, x);
      return p;
    }

    if (type == "op"){
      VectorOrderParamProperty *p = new VectorOrderParamProperty(*d_sys, d_pbc);
      p->parse(arguments, x);
      return p;
    }

    if (type == "pr"){
      PseudoRotationProperty *p = new PseudoRotationProperty(*d_sys, d_pbc);
      p->parse(arguments, x);
      return p;
    }

    if (type == "pa"){
      PuckerAmplitudeProperty *p = new PuckerAmplitudeProperty(*d_sys, d_pbc);
      p->parse(arguments, x);
      return p;
    }

    if (type == "cpq"){
      CPAmplitudeProperty *p = new CPAmplitudeProperty(*d_sys, d_pbc);
      p->parse(arguments, x);
      return p;
    }
    
    if (type == "cpt"){
      CPThetaProperty *p = new CPThetaProperty(*d_sys, d_pbc);
      p->parse(arguments, x);
      return p;
    }
    
    if (type == "cpp"){
      CPPhiProperty *p = new CPPhiProperty(*d_sys, d_pbc);
      p->parse(arguments, x);
      return p;
    }
    
    if (type == "expr"){
      ExpressionProperty *p = new ExpressionProperty(*d_sys, d_pbc);
      p->parse(arguments, x);
      return p;
    }

    if (type == "hb"){
      HBProperty *p = new HBProperty(*d_sys, d_pbc);
      p->parse(arguments, x);
      return p;
    }
    
    if (type == "st"){
      StackingProperty *p = new StackingProperty(*d_sys, d_pbc);
      p->parse(arguments, x);
      return p;
    }

    // throw exception type error
    // either do the user properties first or catch the exception
    std::cerr << "unknown type: " << type << std::endl;
    throw PropertyContainer::Exception(" type unknown\n");
  }

  std::string PropertyContainer::toTitle()const
  {
    std::string s = "";
    // stop output after 10 properties
    const_iterator it = begin();
    for(int c = 0; it != end() && c<10; it++, c++)
      s += (*it)->toTitle() + "\t\t";
    if (it != end())
      s += "[ ... ]\t\t";
      // also print the last one, to acknowledge the range
      s += back()->toTitle();
    return s;
  }
  
  std::string PropertyContainer::toString()const
  {
    std::string s = "";
    for(const_iterator it = begin(); it != end(); it++)
      s += (*it)->toString() + "\t\t";
    return s;
  }
  
  void PropertyContainer::calc()
  {
    for(iterator it = begin(); it != end(); ++it)
      (*it)->calc();
  }
  
  std::ostream &operator<<(std::ostream &os, PropertyContainer const & ps)
  {
    for(PropertyContainer::const_iterator it = ps.begin(); it != ps.end(); it++)
      os << (**it) << "\t\t";
    os << endl;

    return os;
  }
  
}


  











