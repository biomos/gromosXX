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

/**
 * @file VectorSpecifier.cc
 * vector specifier implementation
 */
#include "VectorSpecifier.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <math.h>
#include <sstream>
#include <cstdio>
#include <string>

#include "parse.h"
#include "Value.h"
#include "../bound/Boundary.h"
#include "../gcore/Box.h"

using namespace gcore;

utils::VectorSpecifier::VectorSpecifier(gcore::System &sys, bound::Boundary * pbc, std::string s,
					std::map<std::string, utils::Value> var)
  : d_atomspec(sys),
    d_pbc(pbc)
{
  parse(s, var);
}

utils::VectorSpecifier::VectorSpecifier(utils::VectorSpecifier const & vs)
  : d_atomspec(vs.d_atomspec),
    d_vec(vs.d_vec),
    d_pbc(vs.d_pbc)
{
}


utils::VectorSpecifier::~VectorSpecifier()
{
}

int utils::VectorSpecifier::setSpecifier(std::string s,
					 std::map<std::string, utils::Value> var)
{
  clear();
  parse(s, var);
  return 0;
}

utils::VectorSpecifier & utils::VectorSpecifier::operator=
(
 utils::VectorSpecifier const & vs
 )
{
  d_atomspec = vs.d_atomspec;
  d_vec = vs.d_vec;
  d_pbc = vs.d_pbc;
  
  return *this;
}

gmath::Vec const & utils::VectorSpecifier::operator()()const
{
  if (d_atomspec.size() > 1)
    d_vec = d_atomspec.pos(0) - d_pbc->nearestImage(d_atomspec.pos(0),d_atomspec.pos(1), 
						    d_atomspec.sys()->box());
  else if (d_atomspec.size() == 1)
    d_vec = d_atomspec.pos(0);
  return d_vec;
}

void utils::VectorSpecifier::clear()
{
  d_vec = 0.0;
  d_atomspec.clear();
}

gcore::System & utils::VectorSpecifier::sys()
{
  return *d_atomspec.sys();
}

bound::Boundary * utils::VectorSpecifier::pbc()
{
  return d_pbc;
}

std::string utils::VectorSpecifier::toString()const
{
  std::ostringstream os;
  std::string s = d_atomspec.toString()[0];
  
  if (s!=""){
    os << "atom("  << s << ")";
  }
  else{
    os << "cart(" << d_vec[0]
       << "," << d_vec[1]
       << "," << d_vec[2]
       << ")";
  }
  return os.str();
}

void utils::VectorSpecifier::parse(std::string s,
				   std::map<std::string, utils::Value> & var)
{
  std::string::size_type bra = s.find('(');
  if (bra == std::string::npos)
    throw Exception("wrong format (no ())");
  std::string::size_type ket = find_matching_bracket(s, '(', bra+1);
  if (ket == std::string::npos)
    throw Exception("no closing bracket found!");
  
  std::string b;

  {
    std::istringstream is(s.substr(0, bra));
    is >> b;
  }
  
  std::string rest = s.substr(bra+1, ket - bra - 2);
  
  if (b == "cart"){
    parse_cartesian(rest);
  }
  else if (b == "polar"){
    parse_polar(rest);
  }
  else if (b == "atom"){
    parse_atom(rest, var);
  }
  else{
    throw Exception("wrong format : type '" + b + "'");
  }
}

void utils::VectorSpecifier::parse_atom(std::string s,
					std::map<std::string, utils::Value> & var)
{
  if (var.count("x") > 0){
    // std::cerr << "adding spec with x=" << var["x"] << std::endl;
    d_atomspec.addSpecifier(s, int(rint(var["x"].scalar())));
  }
  else
    d_atomspec.addSpecifier(s);
}

void utils::VectorSpecifier::parse_cartesian(std::string s)
{
  // std::cout << "parse_cartesian : " << s << "..." << std::endl;

  std::string::size_type it = s.find(',');
  if (it == std::string::npos)
    throw Exception("cartesian: vector separated by , expected!");
  
  std::istringstream is(s.substr(0, it));
  if (!(is >> d_vec[0]))
    throw Exception("cartesian: could not read number");
  
  std::string s2 = s.substr(it+1, std::string::npos);

  it = s2.find(',');
  if (it == std::string::npos)
    throw Exception("cartesian: vector separated by , expected!");
  
  is.clear();
  is.str(s2.substr(0, it));
  if (!(is >> d_vec[1]))
    throw Exception("cartesian: could not read number");

  is.clear();
  is.str(s2.substr(it+1, std::string::npos));
  if (!(is >> d_vec[2]))
    throw Exception("cartesian: could not read number");
  
  // std::cout << "cartesian read: " << d_vec[0] << " | " 
  // << d_vec[1] << " | " << d_vec[2] << std::endl;

}


void utils::VectorSpecifier::parse_polar(std::string s)
{
  // std::cout << "parse_polar : " << s << "..." << std::endl;

  d_vec = 0.0;
  double alpha, beta, r;
  
  std::string::size_type it = s.find(',');
  if (it == std::string::npos)
    throw Exception("cartesian: vector separated by , expected!");
  
  std::istringstream is(s.substr(0, it));
  if (!(is >> r))
    throw Exception("cartesian: could not read number");
  
  std::string s2 = s.substr(it+1, std::string::npos);

  it = s2.find(',');
  if (it == std::string::npos)
    throw Exception("cartesian: vector separated by , expected!");
  
  is.clear();
  is.str(s2.substr(0, it));
  if (!(is >> alpha))
    throw Exception("cartesian: could not read number");

  is.clear();
  is.str(s2.substr(it+1, std::string::npos));
  if (!(is >> beta))
    throw Exception("cartesian: could not read number");

  const double cosa = cos(alpha * M_PI / 180.0);
  const double sina = sin(alpha * M_PI / 180.0);
  const double cosb = cos(beta  * M_PI / 180.0);
  const double sinb = sin(beta  * M_PI / 180.0);
  
  d_vec[0] = cosa * cosb * r;
  d_vec[1] = sina * r;
  d_vec[2] = -sinb * cosa * r;

  // std::cout << "vector from polar: " << gmath::v2s(d_vec) << std::endl;

}
