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


#include "Arguments.h"

#include <cstddef>
#include <ostream>
#include <sstream>
#include <fstream>
#include <set>
#include <string>
#include <map>
#include <istream>

#include <config.h>
#include "../utils/debug.h"
#include "../gromos/Exception.h"

using namespace std;

namespace args {

  // Global flags to check if i/o should be in GROMOS96 format
  bool Arguments::inG96 = false;
  bool Arguments::outG96 = false;

  typedef multimap<string, string>::value_type argType;

  // checks if an argument is known

  static int isKnown(const string str, set<string> known_args) {
    // check for g96 mode or develop argument - they are always known
    if (str == "inG96" || str == "outG96" || str == "verb" || str == "develop") return 1;
    
    set<string>::const_iterator fr, to = known_args.end();
    for (fr = known_args.begin(); fr != to; fr++)
      if (*fr == str)break;
    if (fr == to)return 0;
    return 1;
  }

  class Arguments_i {
    friend class Arguments;
    friend istream & operator>>(istream &is, Arguments &args);
    string d_usage;
    string d_prog;
    set<string> d_known;

    Arguments_i() : d_usage(""), d_prog(""), d_known() {
    }

  };

  Arguments::Arguments(int argc, char **argv, const Argument_List & known,
          const string &usage) :
  multimap<string, string>(), d_this(new Arguments_i) {

    if (argc)d_this->d_prog = argv[0];
    d_this->d_usage = "\n#\n" + usage;

    // copy the argument set
    d_this->d_known = known.known;

    string s("");

    for (int i = 1; i < argc; i++) {
      if (string(argv[i]) == "@versioninfo") {
        // remove path from argv[0]
        string program(argv[0]);
        const size_t found = program.find_last_of("/\\");
        program = program.substr(found + 1);

        ostringstream os;
        os << endl
                << "This is GROMOS++ program \"" << program << "\"" << endl
                << "version: " << GROMOS_VERSION << endl
                << "built:   " << GROMOS_DATE << endl;
#ifdef HAVE_GROMOSXX
        os << "GROMOS routines available." << endl;
#endif
        throw gromos::Exception("VERSION INFORMATION", os.str());
      }
    }

    if (argc == 1)
      throw Arguments::Exception(d_this->d_usage);

    for (int i = 1; i < argc; i++) {
      if (string(argv[i]) == "@f") {
        // input file
        ++i;
        ifstream inp(argv[i]);
        if (inp.good() && inp.is_open()) {

          inp >> *this;
          inp.close();
        } else
          throw gromos::Exception("Arguments", "Could not open file " + string(argv[i]));
      } else
        s += string(argv[i]) + ' ';
    }

    // istrstream is(s.c_str());
    stringstream is(s.c_str());
    is >> *this;
    
    // now we have the arguments. Check for g96 mode and remove them from the 
    // argument list.

    { // GROMOS96 input
      iterator it = lower_bound("inG96"), to = upper_bound("inG96");
      if (it != to) { // the argument was found
        inG96 = true;
        erase(it, to);
      }
    }
    { // GROMOS96 output
      iterator it = lower_bound("outG96"), to = upper_bound("outG96");
      if (it != to) { // the argument was found
        outG96 = true;
        erase(it, to);
      }
    } // DEBUG
    {
      iterator it = lower_bound("verb"), to = upper_bound("verb");
      if (it != to) { // do debug
        if (!(istringstream(it->second) >> ::debug_level)) {
          ::debug_level = 0;
          throw gromos::Exception("Arguments", "@verb has to be numeric!");
        }
        erase(it, to);
      }
    }
#ifdef NDEBUG
    if (::debug_level)
      throw gromos::Exception("Arguments", "@verb given but this is a non-debug version of GROMOS++.");
#endif
  }

  Arguments::~Arguments() {
    delete d_this;
  }

  istream & operator>>(istream &istr, Arguments &args) {
    // get away the comments
    string s("");

    while (!istr.eof()) {
      string buff;
      getline(istr, buff);
      s += string(buff);
      if (s.find("#") <= s.length())s = s.erase(s.find("#"));
      s += '\n';
    }
    stringstream is(s.c_str());

    string str, last;

    if (!(is >> last))
      return istr;
    if (last[0] != '@')
      throw Arguments::Exception(args.d_this->d_usage);
    last = last.substr(1);
    if (args.find(last) != args.end())
      args.erase(args.lower_bound(last), args.upper_bound(last));

    if (!isKnown(last, args.d_this->d_known)) {
      string except = "\n#\n# Argument @" + last + " not known! Possible arguments: " + args.d_this->d_usage;
      throw Arguments::Exception(except);
    }

    while (is >> str) {
      if (str == "@help") throw Arguments::Exception(args.d_this->d_usage);

      if (str[0] == '@') {
        if (args.find(last) == args.end())
          args.insert(argType(last, ""));
        last = str.substr(1);
        if (args.find(last) != args.end())
          for (Arguments::iterator l = args.lower_bound(last);
                  l != args.upper_bound(last); ++l)
            args.erase(l);

        if (!isKnown(last, args.d_this->d_known)) {
          string except = "\n\nArgument @" + last + " not known! Possible arguments: " + args.d_this->d_usage;
          throw Arguments::Exception(except);

          continue;
        }
      } else {
        args.insert(argType(last, str));
      }
    }

    // insert the last one without value
    if (args.find(last) == args.end())
      args.insert(argType(last, ""));
    
    return istr;
  }

  const string &Arguments::operator[](const string &str)const {
    const_iterator f = find(str);
    if (f == end()) {
      std::ostringstream os;
      os << "\narguments: could not access '" << str << "'\n";
      os << d_this->d_usage;
      throw Exception(os.str());
    }
    else return find(str)->second;
  }

  int Arguments::check(const string &str, int num_args)const {
    if (find(str) == end())
      throw Exception(d_this->d_usage);
    int num = 0;
    for (const_iterator l = lower_bound(str), u = upper_bound(str);
            l != u; ++l)
      if (l->second != "")++num;
    if (num < num_args)
      throw Exception(d_this->d_usage);
    return 0;
  }

  int Arguments::count(const string &str)const {
    if (find(str) == end())
      return -1;
    int num = 0;
    for (const_iterator l = lower_bound(str), u = upper_bound(str);
            l != u; ++l)
      if (l->second != "")++num;
    return num;
  }

  void Arguments::underDevelopment() {
    if (count("develop") < 0) {
	  throw gromos::Exception(d_this->d_prog, "This program is not tested!\nIf you still want to use it add @develop to the arguments.");
	}
  }

}
