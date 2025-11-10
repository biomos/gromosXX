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
 * @file Value.h
 * Value
 */

/* 	$Id$	 */

#ifndef INCLUDED_UTILS_VALUE
#define INCLUDED_UTILS_VALUE

#include <vector>
#include <string>
#include <sstream>

#include "../gromos/Exception.h"
#include "VectorSpecifier.h"

namespace utils
{
  enum value_enum { val_scalar = 0, val_vector = 1, val_vecspec = 2 };

  class Value
  {
  public:
    Value() 
      : d_value_type(val_scalar), d_scalar(0.0), d_vec(0.0) {}
    explicit Value(double d)
      : d_value_type(val_scalar), d_scalar(d), d_vec(0.0) {}
    explicit Value(gmath::Vec v)
      : d_value_type(val_vector), d_scalar(0.0), d_vec(v) {}
    explicit Value(VectorSpecifier const & v)
      : d_value_type(val_vecspec), d_scalar(0.0), d_vecspec(v) {}
    Value(Value const &v)
      : d_value_type(v.d_value_type), d_scalar(v.d_scalar),
	d_vec(v.d_vec), d_vecspec(v.d_vecspec) {}

    gmath::Vec vec()const
    {
      switch(d_value_type){
	case val_scalar:
	  throw Exception("no vector value");
	case val_vector:
	  return d_vec;
	case val_vecspec:
	  return d_vecspec();
      }
      throw Exception("wrong value type");
    }
    double scalar()const
    {
      if (d_value_type != val_scalar)
	throw Exception("no scalar value");
      return d_scalar;
    }
    
    value_enum type()const
    {
      return d_value_type;
    }

    // maybe they should rather return the Value?
    double operator=(double d)
    {
      d_value_type = val_scalar;
      d_scalar = d;
      return d;
    }
    
    gmath::Vec operator=(gmath::Vec const & v)
    {
      d_value_type = val_vector;
      d_vec = v;
      return d_vec;
    }

    gmath::Vec operator=(VectorSpecifier const & v)
    {
      d_value_type = val_vecspec;
      d_vecspec = v;
      return d_vecspec();
    }
    
    std::string toString()const
    {
      std::ostringstream os;
      switch(d_value_type){
	case val_scalar:
	  os << d_scalar;
	  break;
	case val_vector:
	  os << "(" << d_vec[0] << ", "
	     << d_vec[1] << ", "
	     << d_vec[2] << ")";
	  break;
	case val_vecspec:
	  os << "(" << d_vecspec()[0] << ", "
	     << d_vecspec()[1] << ", "
	     << d_vecspec()[2] << ")";
	  break;
      }
      return os.str();
    }

    void parse(std::string s)
    {
      std::string::size_type bra = s.find('(');
      if (bra != std::string::npos){
	// assume it's a vector specifier
	throw Exception("no system and boundary for VectorSpecifier");
      }
      else{
	std::istringstream is(s);
	if (!(is >> d_scalar))
	  throw Exception("could not parse value");
      }
    }

    void parse(std::string s,
	       std::map<std::string, Value> & var,
	       gcore::System & sys, bound::Boundary *pbc)
    {
      std::string::size_type bra = s.find('(');
      if (bra != std::string::npos){
	// assume it's a vector specifier
	VectorSpecifier vs(sys, pbc, s, var);

	d_value_type = val_vecspec;
	d_vecspec = vs;
      }
      else{
	std::istringstream is(s);
	if (!(is >> d_scalar))
	  throw Exception("could not parse value");
      }
    }

    struct Exception: public gromos::Exception{
      Exception(const std::string &what): 
	gromos::Exception("Value", what){}
    };

    Value & operator-()
    {
      switch(d_value_type){
	case val_scalar:
	  d_scalar = -d_scalar;
	  break;
	case val_vector:
	  d_vec = -d_vec;
	  break;
	case val_vecspec:
	  d_vec = -vec();
	  d_value_type = val_vector;
	  // throw Exception("no unary minus on VectorSpecifier");
      }
      return *this;
    }

    double abs2() const {
      double abs2v = 0.0;
      switch(d_value_type){
	case val_scalar:
          abs2v = d_scalar * d_scalar;
	  break;
	case val_vector:
          abs2v = d_vec.abs2();
	  break;
	case val_vecspec:
          abs2v = vec().abs2();
          break;
      }
      return abs2v;
    }
    
  protected:
    value_enum d_value_type;
    double d_scalar;
    gmath::Vec d_vec;
    VectorSpecifier d_vecspec;
  };
  
  inline std::ostream & operator<<(std::ostream & os, Value const & v)
  {
    os << v.toString();
    return os;
  }

  // some math
  inline Value operator+(Value const &v1, Value const &v2)
  {
    switch(v1.type()){
      case val_scalar:
	{
	  if (v2.type() != val_scalar)
	    throw Value::Exception("adding different types");
	  return Value(v1.scalar() + v2.scalar());
	}
      case val_vector:
      case val_vecspec:
	{
	  if (v2.type() == val_vector ||
	      v2.type() == val_vecspec)
	    return Value(v1.vec() + v2.vec());
	  else
	    throw Value::Exception("adding different types");
	}
    }
    throw Value::Exception("wrong value type");
  }
  inline Value operator-(Value const &v1, Value const &v2)
  {
    switch(v1.type()){
      case val_scalar:
	{
	  if (v2.type() != val_scalar)
	    throw Value::Exception("adding different types");
	  return Value(v1.scalar() - v2.scalar());
	}
      case val_vector:
      case val_vecspec:
	{
	  if (v2.type() == val_vector ||
	      v2.type() == val_vecspec)
	    return Value(v1.vec() - v2.vec());
	  else
	    throw Value::Exception("adding different types");
	}
    }
    throw Value::Exception("wrong value type");
  }
  inline Value operator*(Value const &v1, Value const &v2)
  {
    if (v1.type() == val_scalar){
      if (v2.type() == val_scalar)
	return Value(v1.scalar() * v2.scalar());
      else return Value(v1.scalar() * v2.vec());
    }
    else{
      if (v2.type() == val_scalar)
	return Value(v1.vec() * v2.scalar());
      else throw Value::Exception("use dot / cross to do vector products");
    }
  }
  inline Value operator/(Value const &v1, Value const &v2)
  {
    if (v2.type() != val_scalar)
      throw Value::Exception("division by vector");

    if (v1.type() == val_scalar)
      return Value(v1.scalar() / v2.scalar());
    else
      return Value(v1.vec() / v2.scalar());
  }

  inline Value operator!(Value const & v)
  {
    if (v.type() != val_scalar)
      throw Value::Exception("operator not (!) not defined for non-scalar values");
    return Value(!v.scalar());
  }

  inline Value operator==(Value const & v1, Value const & v2)
  {
    if (v1.type() != v2.type())
      throw Value::Exception("comparing values of unequal type not allowed");
    
    if (v1.type() == val_scalar)
      return Value(v1.scalar() == v2.scalar());
    
    return Value(v1.vec() == v2.vec());
  }

  inline bool operator==(Value const & v, bool b)
  {
    if (v.type() != val_scalar)
      throw Value::Exception("comparison with bool only allowed for scalar values");
    
    return (v.scalar() != 0);
  }

  inline Value operator!=(Value const & v1, Value const & v2)
  {
    if (v1.type() != v2.type())
      throw Value::Exception("comparing values of unequal type not allowed");
    
    if (v1.type() == val_scalar)
      return Value(v1.scalar() != v2.scalar());
    
    return !Value(v1.vec() == v2.vec());
  }

  inline Value operator<(Value const & v1, Value const & v2)
  {
    if (v1.type() != v2.type())
      throw Value::Exception("comparing values of unequal type not allowed");
    
    if (v1.type() == val_scalar)
      return Value(v1.scalar() < v2.scalar());
    
    return Value(v1.vec().abs2() < v2.vec().abs2());
  }

  inline Value operator>(Value const & v1, Value const & v2)
  {
    if (v1.type() != v2.type())
      throw Value::Exception("comparing values of unequal type not allowed");
    
    if (v1.type() == val_scalar)
      return Value(v1.scalar() > v2.scalar());
    
    return Value(v1.vec().abs2() > v2.vec().abs2());
  }

  inline Value operator<=(Value const & v1, Value const & v2)
  {
    if (v1.type() != v2.type())
      throw Value::Exception("comparing values of unequal type not allowed");
    
    if (v1.type() == val_scalar)
      return Value(v1.scalar() <= v2.scalar());
    
    return Value(v1.vec().abs2() <= v2.vec().abs2());
  }

  inline Value operator>=(Value const & v1, Value const & v2)
  {
    if (v1.type() != v2.type())
      throw Value::Exception("comparing values of unequal type not allowed");
    
    if (v1.type() == val_scalar)
      return Value(v1.scalar() >= v2.scalar());
    
    return Value(v1.vec().abs2() >= v2.vec().abs2());
  }

  inline Value operator&&(Value const & v1, Value const & v2)
  {
    if (v1.type() != v2.type())
      throw Value::Exception("comparing values of unequal type not allowed");
    
    if (v1.type() == val_scalar)
      return Value(v1.scalar() && v2.scalar());

    throw Value::Exception("operator AND (&&) of non-scalars not allowed");
  }

  inline Value operator||(Value const & v1, Value const & v2)
  {
    if (v1.type() != v2.type())
      throw Value::Exception("comparing values of unequal type not allowed");
    
    if (v1.type() == val_scalar)
      return Value(v1.scalar() || v2.scalar());

    throw Value::Exception("operator OR (||) of non-scalars not allowed");
  }

  inline Value abs(Value const &v)
  {
    if (v.type() == val_scalar)
      return Value(fabs(v.scalar()));
    return Value(v.vec().abs());
  }
  inline Value abs2(Value const &v)
  {
    if (v.type() == val_scalar)
      return Value(v.scalar() * v.scalar());
    return Value(v.vec().abs2());
  }
  inline Value dot(Value const &v1, Value const &v2)
  {
    if (v1.type() == val_scalar ||
	v2.type() == val_scalar)
      throw Value::Exception("dot product requires two vectors as arguments");

    return Value(v1.vec().dot(v2.vec()));
  }
  inline Value cross(Value const &v1, Value const &v2)
  {
    if (v1.type() == val_scalar ||
	v2.type() == val_scalar)
      throw Value::Exception("cross product requires two vectors as arguments");

    return Value(v1.vec().cross(v2.vec()));
  }
  inline Value sin(Value const &v1)
  {
    if (v1.type() == val_scalar)
      return Value(::sin(v1.scalar()));
    throw Value::Exception("sin() only defined for scalar argument");
  }
  inline Value cos(Value const &v1)
  {
    if (v1.type() == val_scalar)
      return Value(::cos(v1.scalar()));
    throw Value::Exception("cos() only defined for scalar argument");
  }
  inline Value tan(Value const &v1)
  {
    if (v1.type() == val_scalar)
      return Value(::tan(v1.scalar()));
    throw Value::Exception("tan() only defined for scalar argument");
  }
  inline Value asin(Value const &v1)
  {
    if (v1.type() == val_scalar)
      return Value(::asin(v1.scalar()));
    throw Value::Exception("asin() only defined for scalar argument");
  }
  inline Value acos(Value const &v1)
  {
    if (v1.type() == val_scalar)
      return Value(::acos(v1.scalar()));
    throw Value::Exception("acos() only defined for scalar argument");
  }
  inline Value atan(Value const &v1)
  {
    if (v1.type() == val_scalar)
      return Value(::atan(v1.scalar()));
    throw Value::Exception("atan() only defined for scalar argument");
  }
  inline Value log(Value const &v1)
  {
    if (v1.type() == val_scalar)
      return Value(::log(v1.scalar()));
    throw Value::Exception("log() only defined for scalar argument");
  }
  inline Value sqrt(Value const &v1)
  {
    if (v1.type() == val_scalar)
      return Value(::sqrt(v1.scalar()));
    throw Value::Exception("sqrt() only defined for scalar argument");
  }
  inline Value exp(Value const &v1)
  {
    if (v1.type() == val_scalar)
      return Value(::exp(v1.scalar()));
    throw Value::Exception("exp() only defined for scalar argument");
  }

  inline long double sin(long double d)
  {
    return ::sin(d);
  }
  inline long double cos(long double d)
  {
    return ::cos(d);
  }
  inline long double tan(long double d)
  {
    return ::tan(d);
  }
  inline long double asin(long double d)
  {
    return ::asin(d);
  }
  inline long double acos(long double d)
  {
    return ::acos(d);
  }
  inline long double atan(long double d)
  {
    return ::atan(d);
  }
  inline long double sqrt(long double d)
  {
    return ::sqrt(d);
  }
  inline long double exp(long double d)
  {
    return ::exp(d);
  }
  inline long double log(long double d)
  {
    return ::log(d);
  }
  inline int abs(int val)
  {
    return std::abs(val);
  }
  inline double abs(double d)
  {
    return ::fabs(d);
  }
  inline long double abs(long double d)
  {
    return std::abs(d);
  }
  
}

#endif

