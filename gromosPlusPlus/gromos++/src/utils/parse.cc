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
#include "parse.h"

#include <map>
#include <string>
#include <iostream>
#include <cassert>
#include <vector>

#include "ExpressionParser.h"
#include "../gromos/Exception.h"

namespace utils
{
  
  std::string::size_type find_par
  (
   std::string s,
   char c,
   std::string::size_type it,
   std::string bra,
   std::string ket
   )
  {
    // recognized brackets: (, [, {, <
    int level = 0;
    
    for( ; it < s.length(); ++it){
      
      if (bra.find(s[it]) != std::string::npos) ++level;
      if (ket.find(s[it]) != std::string::npos) --level;
      
      
      if (s[it] == c && level == 0) return it;
      
    }
    
    return std::string::npos;
  }

  std::string::size_type find_par
  (
   std::string s,
   string c,
   std::string::size_type it,
   std::string bra,
   std::string ket
   )
  {
    // recognized brackets: (, [, {, <
    int level = 0;
    
    for( ; it < s.length(); ++it){
      
      if (bra.find(s[it]) != std::string::npos) ++level;
      if (ket.find(s[it]) != std::string::npos) --level;

      
      if (c.find(s[it]) != std::string::npos && level == 0) return it;
    }
    
    return std::string::npos;
  }
  
  std::string::size_type find_matching_bracket
  (
   std::string s,
   char bra,
   std::string::size_type it
   )
  {
    char ket;
    if (bra == '(') ket = ')';
    else if (bra == '[') ket = ']';
    else if (bra == '{') ket = '}';
    else if (bra == '<') ket = '>';
    else{
      std::cerr << "could not determine brackets" << std::endl;
      throw gromos::Exception("parser", "Bracket not recognised");
    }
    
    int level = 1;
    for( ; it < s.length() && level != 0; ++it){
      if (s[it] == bra) ++level;
      else if (s[it] == ket) --level;
    }
    
    if (level) return std::string::npos;
    else return it;
  }

  /**
   * parse a range into an
   * index array
   */
  void parse_range(std::string s, std::vector<int> & range, int x)
  {
    std::map<std::string, int> var;
    var["x"] = x;
    
    ExpressionParser<int> ep;
    
    std::string::size_type it = s.find(',');

    std::string rest;
    if (it != std::string::npos){
      rest = s.substr(it+1, std::string::npos);
      s = s.substr(0, it);
    }
    else{
      rest = "";
    }

    // parse the single part
    it = s.find('-');
    //  +-- Check for this because of ranges like "-1"
    //  |
    if (it && it != std::string::npos){
      int sr, er;
      sr = ep.parse_expression(s.substr(0, it), var);
      er = ep.parse_expression(s.substr(it+1, std::string::npos), var);
      for(int i=sr; i<=er; ++i)
	range.push_back(i);
    }
    else{
      int i = ep.parse_expression(s, var);
      range.push_back(i);
    }

    if (rest != "")
      parse_range(rest, range, x);
  }
}
