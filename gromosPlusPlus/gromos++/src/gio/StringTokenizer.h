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

// StringTokenizer.h

#ifndef INCLUDED_GIO_STRINGTOKENIZER
#define INCLUDED_GIO_STRINGTOKENIZER

#include <string>
#include <vector>

#include "../gromos/Exception.h"


namespace gio{
  
  /**
   * Class StringTokenizer
   * Takes a string and returns tokens as vector of strings
   * separated by a delimiter.
   *
   * @class StringTokenizer
   * @ingroup gio
   * @author M.A. Kastenholz
   */
  class StringTokenizer{
    // not implemented
    StringTokenizer(const StringTokenizer&);
    StringTokenizer &operator=(const StringTokenizer&);
    std::string to_tokenize;
    std::string delimiter;

    public:
    // Constructors
    StringTokenizer(const std::string &str, const std::string& delim = " ");
    ~StringTokenizer() {};

    // Methods
    void setdelimiter(const std::string delim);
    void setstring(const std::string str);
    std::vector<std::string> tokenize();
    
    //Exceptions
    struct Exception: public gromos::Exception{
      Exception(const std::string& what_arg) : gromos::Exception("StringTokenizer", what_arg){}
    };
  };
}
#endif
