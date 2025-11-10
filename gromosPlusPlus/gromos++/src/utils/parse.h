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


#ifndef INCLUDED_UTILS_PARSE
#define INCLUDED_UTILS_PARSE

#include <string>
#include <vector>

namespace utils
{

  /**
   * find matching bracket
   */
  std::string::size_type find_matching_bracket
  (
   std::string s,
   char bra='(',
   std::string::size_type it=0
   );

  /**
   * find character, taking care of brackets in string
   * ie: looking for ';' in va(prop1;prop2);prop3
   * will find second ';'
   */
  std::string::size_type find_par
  (
   std::string s,
   char c=';',
   std::string::size_type it=0,
   std::string bra = "([{<",
   std::string ket = ">}])"
   );

  /**
   * find character contained in c, taking care of brackets in string
   * ie: looking for ';' in va(prop1;prop2);prop3
   * will find second ';'
   */
  std::string::size_type find_par
  (
   std::string s,
   std::string c,
   std::string::size_type it=0,
   std::string bra = "([{<",
   std::string ket = ">}])"
   );

  /**
   * parse a range into an
   * index array
   */
  void parse_range(std::string s, std::vector<int> & range, int x=-1);
  
}


#endif
