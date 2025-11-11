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

// gromos_exception.h

#ifndef INCLUDED_GROMOS_EXCEPTION
#define INCLUDED_GROMOS_EXCEPTION

#include <exception>
#include <string>

namespace gromos
{
  /**
   * Class Exception
   * Gromos++ exception class.
   *
   * Recently the use of this class has been questioned. Now that I write 
   * the description I sort of seem to agree.
   *
   * @class Exception
   * @author R. Buergi
   * @ingroup gromos
   */
  class Exception: public std::exception{
    std::string d_mesg;
  public:
    /**
     * Exception constructor
     * @param from A description of which program / routine 
     *             throws the description
     * @param mess The error message
     */
    Exception(const std::string from,
	      const std::string mess) throw():
      d_mesg(from + ": " + mess) 
      {}
  virtual  const char *what() const throw(){
      return d_mesg.c_str();
    }
  virtual  ~Exception() throw(){}
  };
}



#endif
