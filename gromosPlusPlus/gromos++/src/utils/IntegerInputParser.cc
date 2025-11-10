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
#include "IntegerInputParser.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include "../gromos/Exception.h"

void utils::IntegerInputParser::parse(std::string const s, int maxnum)
{
  if(s=="ALL" || s=="all"){
    for(int i=0; i<maxnum; i++) insert(i+1);
    return;
  }
  std::string::size_type iterator;
  if((iterator=s.find(',')) != std::string::npos){
    parse(s.substr(0,iterator), maxnum);
    parse(s.substr(iterator+1, std::string::npos), maxnum);
  }
  else{
    std::istringstream is;
    int rangeBegin, rangeEnd;
    if((iterator=s.find('-')) != std::string::npos){
      is.str(s.substr(0,iterator));
      if(!(is >> rangeBegin))
	throw gromos::Exception("IntegerInputParser", 
				"Invalid begin of range "+ s);
      is.clear();
      is.str(s.substr(iterator+1, std::string::npos));
      if(!(is >> rangeEnd))
	throw gromos::Exception("IntegerInputParser", 
				"Invalid end of range " + s);
      for(int i=rangeBegin; i<= rangeEnd; ++i){
	if(i> maxnum)
	  throw gromos::Exception("IntegerInputParser",
				  "Requested number too high: "+s);
	insert(i);
      }
    }
    else{
      is.str(s);
      if(!(is >> rangeBegin))
	throw gromos::Exception("IntegerInputParser", 
				"Invalid number specified "+ s);
      if(rangeBegin > maxnum)
	throw gromos::Exception("IntegerInputParser",
				"Requested number too high: "+s);
      insert(rangeBegin);
    }
  }
}

void utils::IntegerInputParser::addSpecifier(std::string const s, int maxnum)
{
  parse(s, maxnum);
}
