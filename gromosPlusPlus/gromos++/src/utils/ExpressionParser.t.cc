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
#include "ExpressionParser.h"

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <sstream>
#include <cassert>
#include <cstdlib>

#include "Value.h"
#include "../gio/InTopology.h"
#include "../gcore/System.h"
#include "../gio/InG96.h"
#include "../bound/Boundary.h"
#include "../args/Arguments.h"
#include "../args/BoundaryParser.h"
#include "../args/GatherParser.h"
#include "../gromos/Exception.h"

int debug_level = 0;

int main(int argc, char *argv[]) {
  std::string usage = argv[0];
  usage += "\n";
  usage += "\t@topo   <topology>\n";
  usage += "\t@pbc    <pbc>\n";
  usage += "\t@traj   <trajectory>\n";
  usage += "\t@type   <expression type>\n";
  usage += "\t@expr   <expression>\n";
  usage += "\t@expr2  <expression>\n";
  usage += "\t@expr3  <expression>\n";
  usage += "\t@expr4  <expression>\n";

  args::Argument_List knowns;
  knowns << "topo" << "pbc" << "traj" << "type" << "expr" << "expr2" << "expr3" << "expr4";
  args::Arguments args(argc, argv, knowns, usage);

  try {

    std::vector<std::string> exprstr;
    {
      std::string s = "";
      for (args::Arguments::const_iterator
        iter = args.lower_bound("expr"),
              to = args.upper_bound("expr");
              iter != to; ++iter) {
        s += " " + iter->second;
      }
      exprstr.push_back(s);
    }


    if (args.count("expr2") > 0) {
      std::string s = "";
      for (args::Arguments::const_iterator
        iter = args.lower_bound("expr2"),
              to = args.upper_bound("expr2");
              iter != to; ++iter) {
        s += " " + iter->second;
      }
      exprstr.push_back(s);
    }
    if (args.count("expr3") > 0) {
      std::string s = "";
      for (args::Arguments::const_iterator
        iter = args.lower_bound("expr3"),
              to = args.upper_bound("expr3");
              iter != to; ++iter) {
        s += " " + iter->second;
      }
      exprstr.push_back(s);
    }
    if (args.count("expr4") > 0) {
      std::string s;
      for (args::Arguments::const_iterator
        iter = args.lower_bound("expr4"),
              to = args.upper_bound("expr4");
              iter != to; ++iter) {
        s += " " + iter->second;
      }
      exprstr.push_back(s);
    }

    if (exprstr.size() < 1)
      throw args::Arguments::Exception("no expressions specified!");

    if (args["type"] == "double") {
      std::cout << "parsing 'double' expression" << std::endl;

      utils::ExpressionParser<double> ep;

      std::map<std::string, double> var;
      var["bla"] = 0.12345;

      try {
        double d = ep.parse_expression(exprstr[0], var);

        std::cout << "result = " << d << std::endl;
      }      catch (gromos::Exception e) {
        std::cerr << "runtime error: " << e.what() << std::endl;
      }
    }
    else if (args["type"] == "int") {
      std::cout << "parsing 'int' expression\n";
      // std::cout << "\t" << exprstr << std::endl;

      utils::ExpressionParser<int> ep;

      std::map<std::string, int> var;
      var["bla"] = 12345;

      std::map<std::string, std::vector<utils::ExpressionParser<int>::expr_struct> > expr_map;

      try {
        for (unsigned int i = 0; i < exprstr.size(); ++i) {

          std::vector<utils::ExpressionParser<int>::expr_struct> expr;

          std::ostringstream os;
          os << "e" << i;

          ep.parse_expression(exprstr[i], var, expr);
          expr_map[os.str()] = expr;

          ep.print_expression(expr);
        }

        ep.calculate(expr_map, var);

        for (unsigned int i = 0; i < exprstr.size(); ++i) {
          std::ostringstream os;
          os << "e" << i;
          std::cout << os.str() << " = " << var[os.str()] << std::endl;
        }
      }      catch (gromos::Exception e) {
        std::cerr << "runtime error: " << e.what() << std::endl;
      }      catch (std::string n) {
        std::cerr << "variable '" + n + "' unknown" << std::endl;
      }

    } else {

      std::cout << "parsing 'value' expression\n";

      try {
        gio::InTopology it(args["topo"]);
        gcore::System sys(it.system());
        gcore::System refSys(it.system());
        bound::Boundary *pbc = args::BoundaryParser::boundary(sys, args);

        // parse gather method
        bound::Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, refSys, args);

        utils::ExpressionParser<utils::Value> ep(&sys, pbc);

        std::map<std::string, utils::Value> var;
        var["bla"] = utils::Value(0.12345);

        gio::InG96 ic(args["traj"]);

        // ic.select("ALL");
        ic >> sys;
        (*pbc.*gathmethod)();

        ic.close();

        try {
          std::cerr << "trying to parse " << exprstr[0] << std::endl;
          utils::Value v = ep.parse_expression(exprstr[0], var);

          std::cout << "result = " << v << std::endl;
        } catch (gromos::Exception e) {
          std::cerr << "runtime error: " << e.what() << std::endl;
        }
      } catch (gromos::Exception const & e) {
        std::cerr << "Exception:\t";
        std::cerr << e.what() << std::endl;
        return 1;
      }
    }
  } catch (gromos::Exception const & e) {
    std::cerr << "Exception:\t";
    std::cerr << e.what() << std::endl;
    return 1;
  }

  return 0;
}
