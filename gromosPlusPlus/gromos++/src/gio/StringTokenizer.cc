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

// StringTokenizer.cc
#include "StringTokenizer.h"

#include <string>
#include <vector>

using gio::StringTokenizer;

  StringTokenizer::StringTokenizer(const std::string &str, const std::string& delim) {
    to_tokenize = str;
    delimiter = delim;
  };
  


  void StringTokenizer::setdelimiter(const std::string delim) {
   delimiter = delim;
  }
  void StringTokenizer::setstring(const std::string str) {
   to_tokenize = str;
  }

  std::vector<std::string> StringTokenizer::tokenize() {
     std::vector< std::string> tokens;
   // Skip delimiters at beginning.
  std::string::size_type lastPos = to_tokenize.find_first_not_of(delimiter, 0);
  // Find first "non-delimiter".
  std::string::size_type pos     = to_tokenize.find_first_of(delimiter, lastPos);

  while (std::string::npos != pos || std::string::npos != lastPos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(to_tokenize.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = to_tokenize.find_first_not_of(delimiter, pos);
      // Find next "non-delimiter"
      pos = to_tokenize.find_first_of(delimiter, lastPos);
    }

   return tokens;
 
  }


