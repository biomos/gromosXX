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
 * @file CommandLine.h
 * Reading stuff from the command line in a fancy way
 */

#ifndef INCLUDED_UTILS_COMMANDLINE
#define INCLUDED_UTILS_COMMANDLINE

#include <iostream>
#include <string>
#include <sstream>
#include <set>
#include <cstdio>

namespace utils {
  namespace CommandLine {
    /**
     * read a line using cin or the readline interface
     * @param[in] prompt the prompt
     * @param[in] out the ostream to write the prompt to.
     * @return the line
     */
    std::string getLine(const std::string & prompt = "",
            std::ostream & out = std::cerr);

    /**
     * Get a value of type T from the command line
     * @param[out] val the value read
     * @param[in] prompt the prompt
     * @param[in] out the ostream to write the prompt to.
     * @param[in] error the error message issues on failure
     */
    template<typename T>
    void getValue(T & val,
            const std::string & prompt = "",
            std::ostream & out = std::cerr,
            const std::string & error = "This is an invalid value. Try again.\n") {
      bool ok = false;
      do {
        std::istringstream s(getLine(prompt, out));
        s >> val;
        if (s.fail())
          out << error;
        else
          ok = true;
      } while (!ok);
    }

    /**
     * Get a value a yes or no answer
     * @param[in] prompt the prompt
     * @param[in] out the ostream to write the prompt to.
     * @param[in] error the error message issues on failure
     * @return true if yes, false if no
     */
    bool getYesNo(const std::string & prompt = "",
            std::ostream & out = std::cerr,
            const std::string & error = "This is an invalid value. Give y(es) or n(o).\n");

    /**
     * Get a value from a given set
     * @param[out] val the value read
     * @param[in] set the set
     * @param[in] prompt the prompt
     * @param[in] out the ostream to write the prompt to.
     * @param[in] error the error message issues on failure. %s is replaced by the set values.
     */
    template<typename T>
    void getValueFromSet(T & val,
            const std::set<T> & set,
            const std::string & prompt = "",
            std::ostream & out = std::cerr,
            const std::string & error = "This is an invalid value.\nValid values are %s\nTry again.\n") {
      std::ostringstream setvalues;
      for(typename std::set<T>::const_iterator it = set.begin(), to = set.end(); it != to; ++it)
        setvalues << *it << " ";

      char buf[error.length() + setvalues.str().length()];
      sprintf(buf, error.c_str(), setvalues.str().c_str());
      std::string error_msg(buf);

      bool ok = false;
      do {
        std::istringstream s(getLine(prompt, out));
        s >> val;
        if (s.fail() || set.find(val) == set.end())
          out << error_msg;
        else
          ok = true;
      } while (!ok);
    }
    
  } // ns commandline
} // ns utils

#endif	/* INCLUDED_UTILS_COMMANDLINE */

