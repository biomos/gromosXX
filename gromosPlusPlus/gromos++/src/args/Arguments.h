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

#ifndef INCLUDED_ARGS_ARGUMENTS
#define INCLUDED_ARGS_ARGUMENTS

#include <map>
#include <set>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <typeinfo>

#include "../gromos/Exception.h"


namespace args {
  class Arguments_i;

  /**
   * Class Argument_List
   * Purpose: add an easy addition interface for the arguments set
   *
   * @class Argument_List
   * @version $Date: Mon Sep 16 14:27:41 MEST 2008
   * @author  N. Schmid
   * @ingroup args
   * @sa args::Arguments
   */
  class Argument_List {
  public:
    std::set<std::string> known;

    /**
     * adds an argument (c string) to the argument list
     */
    Argument_List & operator<<(const char* arg) {
      known.insert(std::string(arg));
      return *this;
    }
  };

    /**
   * Class Arguments
   * Purpose: Parse arguments from the command line or an input file.
   *
   * It also takes care of GROMOS96 compatibility: The arguments \@inG96 and
   * \@outG96 are known to all programs.
   *
   *
   * Description:
   * This class is used to parse arguments from the command line or an input file .
   *
   *
   * @class Arguments
   * @version $Date: Mon Jul 15 14:17:41 MEST 2002
   * @author  R. Buergi
   * @author  M. Kastenholz
   * @ingroup args
   * @sa args::BoundaryParser
   * @sa args::GatherParser
   */

  class Arguments : public std::multimap<std::string, std::string> {
    Arguments_i *d_this;
    // not implemented
    Arguments();
    Arguments(const Arguments &);
    Arguments & operator=(const Arguments &);
  public:
    /**
     * this variable determines whether GROMOS++ is in GROMOS96 input format
     * mode
     */
    static bool inG96;
    /**
     * this variable determines whether GROMOS++ is in GROMOS96 output format
     * mode
     */
    static bool outG96;
    /**
     * Arguments constructor.
     * Details.
     */
    Arguments(int argc, char **argv, const Argument_List & known, const std::string &usage);
    /**
     * Arguments destructor.
     * Details.
     */
    ~Arguments();

    /**
     * Checks for whether the argument string had num_arg arguments
     * @param &str Takes a std::string as argument.
     * @param num_args Integer of the number of arguments.
     * @return check Integer to check for failure.
     */
    int check(const std::string &str, int num_args = 0)const;

    /**
     * Returns the number of arguments that follow string
     * @param &str Takes a std::string as argument.
     * @return the number of arguments found for this argument. Returns -1
     * if string was not found at all in the argument list.
     */
    int count(const std::string &str)const;

    // This has to be in to fix a bug in gcc Solaris 2.6 (?)
    //   const const_iterator begin()const
    //   {return this->multimap<string,string>::begin();}

    /**
     * Member operator >>.
     * Details.
     */
    friend std::istream & operator>>(std::istream &is, Arguments &args);

    /**
     * Member operator [] used to access the members.
     * Details.
     */
    const std::string & operator[](const std::string &str)const;

    class Defaults;

    /**
     * Get a value, like a double or an int from an argument
     * @param arg the name of the argument
     * @param required is the argument required or is there a default
     * @param def the default value if no argument is given
     */
    template<typename T>
    T getValue(const std::string & arg, bool required = true, const T & def = 0);

    /**
     * Get multiple values as a vector
     * @param arg the name of the argument
     * @param num the number of arguments
     * @param required is the argument required or is there a default
     * @param def a vector containing the default values
     */
    template<typename T>
    std::vector<T> getValues(const std::string & arg, unsigned int num, bool required = true, const std::vector<T> & def = std::vector<T>());

        /**
     * @struct Exception
     * Throws the Usage string if invoked.
     */
    struct Exception : public gromos::Exception {
      /**
       * @exception const std::string &str throws Usage string if invoked.
       */
      Exception(const std::string & str) : gromos::Exception("# Usage", str) {
      }
    };

    /**
     * a little helper class for argument defaults
     * @class Default
     */
    template<typename T>
    class Default : public std::vector<T> {
    public:
      Default & operator<<(const T & val) {
        this->push_back(val); return *this;
      }
    };

	/**
	* Checks if the program is under development and crashes if so (unless
	* there is the argument \@develop)
	*/
    void underDevelopment();

  };

  template<typename T>
  T Arguments::getValue(const std::string & arg, bool required, const T & def) {
    Arguments::const_iterator it = lower_bound(arg), to = upper_bound(arg);
    if (it != to) {
      std::istringstream is(it->second);
      T val;
      if (is >> val) {
        return val;
      } else {
        std::ostringstream msg;
        msg << "Argument '" << arg << "' has an invalid value." << std::endl
                << "#        Cannot interpret '" << it->second << "' as type "
                << typeid(T).name() << "." << std::endl;
        throw Arguments::Exception(msg.str());
      }
      if (++it != to) {
        std::ostringstream msg;
        msg << "Too many values for argument '" << arg << "'.";
        throw Arguments::Exception(msg.str());
      }
    }
    if (required) {
      std::ostringstream msg;
      msg << "Argument '" << arg << "' is required.";
      throw Arguments::Exception(msg.str());
    }

    return def;
  }

  template<typename T>
  std::vector<T> Arguments::getValues(const std::string & arg, unsigned int num, bool required, const std::vector<T> & def) {
    Arguments::const_iterator it = lower_bound(arg), to = upper_bound(arg);

    if (it == to) {
      if (required) {
        std::ostringstream msg;
        msg << "Argument '" << arg << "' is required (" << num << " values).";
        throw Arguments::Exception(msg.str());
      }
      if (def.size() != num) {
        throw Arguments::Exception("Not enough default arguments provided.");
      }
      return def;
    }
    std::vector<T> result(num);
    for (unsigned int i = 0; i < num; ++i, ++it) {
      if (it == to) {
        std::ostringstream msg;
        msg << "Not enough values for argument '" << arg << "' given.";
        throw Arguments::Exception(msg.str());
      }
      std::istringstream is(it->second);
      if (!(is >> result[i])) {
        std::ostringstream msg;
        msg << "Argument '" << arg << "' has an invalid value (value number "
                << i+1 << ")." << std::endl
                << "#        Cannot interpret '" << it->second << "' as type "
                << typeid(T).name() << "." << std::endl;
        throw Arguments::Exception(msg.str());
      }
    }
    if (it != to) {
      std::ostringstream msg;
      msg << "Too many arguments for argument '" << arg << "'.";
      throw Arguments::Exception(msg.str());
    } else {
      return result;
    }
  }
}

#endif







