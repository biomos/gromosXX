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
#include "OutformatParser.h"

#include <cctype>
#include <string>
#include <sstream>
#include <cassert>
#include <iostream>
#include <algorithm>

#include "Arguments.h"
#include "../gromos/Exception.h"
#include "../gio/OutCoordinates.h"
#include "../gio/OutG96S.h"
#include "../gio/OutG96.h"
#include "../gio/OutPdb.h"
#include "../gio/Outvmdam.h"
#include "../gio/OutCif.h"

using namespace std;
using namespace gcore;
using namespace args;
using namespace gio;

OutCoordinates * args::OutformatParser::parse(Arguments & args,
        string & ext, string argname) {
  OutCoordinates *oc;
  ext = ".cnf";
  if (args.count("outformat") > 0) {
      Arguments::const_iterator it = args.lower_bound("outformat"),
                                to = args.upper_bound("outformat");
      string format = args["outformat"];
      transform( format.begin(), format.end(),
                 format.begin(),
                 [] (unsigned char c) { return std::tolower(c); } );
      if ((format == "pdb") || (format == "pqr")){
          ++it;

          std::string flavour = format;

          if (it == to) {
              oc = new OutPdb(flavour);
          } else {
              double factor = 10.0;
              bool renumber=false;
              while (it != to) {
                  istringstream is(it->second);
                  if (is.str() == "renumber") { 
                      renumber=true;
                  }
                  else if (!(is >> factor))
                      throw gromos::Exception("OutformatParser", "@outformat pdb factor has to be numeric.!");
                  ++it;
              }
              oc = new OutPdb(flavour, factor, renumber);
          }
          ext = "." + flavour;  //.pdb or .pqr
      } else if (format == "cif") {
          ++it;

          if (it == to) {
              oc = new OutCif();
          } else {
              double factor = 10.0;
              bool renumber=false;
              while ( it != to ) {
                  istringstream is( it->second );
                  if (is.str() == "renumber") { 
                      renumber=true;
                  } else if (!(is >> factor)) {
                      throw gromos::Exception("OutformatParser", "@outformat mmCIF factor has to be numeric.!");
                  }
                  ++it;
              }
              oc = new OutCif( factor, renumber );
          }
          ext = ".cif";
      } else if (format == "cnf") {
          oc = new OutG96S();
          ext = ".cnf";
      } else if (format == "trc") {
          oc = new OutG96();
          ext = ".trc";
      } else if (format == "por") {
          oc = new OutG96S(true);
          ext = ".por";
      } else if (format == "vmdam") {
          ++it;
          if (it == to) {
              oc = new Outvmdam();
          } else {
              istringstream is(it->second);
              double factor = 10.0;
              if (!(is >> factor))
                  throw gromos::Exception("OutformatParser", "@outformat vmdam factor has to be numeric.!");
              oc = new Outvmdam(factor);
          }
          ext = ".vmd";
      } else {
          ostringstream msg;
          msg << "Output format '" << format << "' is unkown." << endl
              << "Known formats are:" << endl
              << "    - cnf" << endl
              << "      Configuration format containing the POSITION block." << endl
              << "    - trc" << endl
              << "      Trajectory format containing the POSITIONRED block." << endl
              << "    - por" << endl
              << "      Position restraints specification format." << endl
              << "    - pdb [<factor to convert length unit to Angstrom, 10.0>]" << endl
              << "      Protein Data Bank (PDB) format." << endl
              << "    - pqr [<factor to convert length unit to Angstrom, 10.0>]" << endl
              << "      Modified Protein Data Bank (PDB) format." << endl
              << "    - vmdam [<factor to convert length unit to Angstrom, 10.0>]" << endl
              << "      VMD's Amber Coordinates format." << endl
              << "    - cif [<factor to convert length unit to Angstrom, 10.0>]" << endl
              << "      crystallographic information file (mmCIF) format." << endl;

          throw gromos::Exception("OutformatParser", msg.str());
      }
  } else {
      oc = new OutG96S();
      ext = ".cnf";
  }

  return oc;
}


