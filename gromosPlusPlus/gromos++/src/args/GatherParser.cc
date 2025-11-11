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

// args_GatherParser.cc

#include "GatherParser.h"

#include <cassert>
#include <iostream>
#include <ostream>
#include <string>

#include "BoundaryParser.h"
#include "Arguments.h"
#include "../gcore/System.h"
#include "../bound/Boundary.h"
#include "../utils/AtomSpecifier.h"
#include "../gromos/Exception.h"

using namespace args;
using namespace bound;
using namespace utils;
using namespace gcore;

bound::Boundary::MemPtr GatherParser::parse(gcore::System &sys, gcore::System &refSys,
					    const Arguments &gathargs,
					    const std::string &str,
					    const bool writeout) {

  string usage = "\n\nAvailable gathering methods\n";
  usage += "\n\t1 (or glist) : based on a list (default). If no list is given, a warning message will be shown";
  usage += "\n\t2 (or gtime) : based on previous frame (time)";
  usage += "\n\t3 (or gref)  : based on a reference structure (refg required)";
  usage += "\n\t4 (or gltime): the 1st frame based on a list, then previous frame (list+time)";
  usage += "\n\t5 (or grtime): the 1st frame based on a ref, then previous frame (ref+time, refg required)";
  usage += "\n\t6 (or gbond) : based on bond connectivity (useful for a single molecule)";
  usage += "\n\t7 (or cog) : based on the center of geometry of the first molecule";
  usage += "\n\t8 (or gfit) : in the 1st frame gather selected molecules to a superimposed reference, in subsequent frames gather based on previous frame; non-selected molecules are gathered to cog of selected molecules";
  usage += "\n\t0 (or nog)   : no gathering\n";

  usage += "\nThree sub-options are available:";
  usage += "\n\t list   [atom_list (optional for methods 1 and 4)]";
  usage += "\n\t refg   <reference molecule for gathering (required for methods 3 and 5)>";
  usage += "\n\t molecules   <optional with method gfit: specify molecules which will be gathered to the reference;\n all other molecules will be gathered with respect to the cog of the selected molecules>\n";

  Boundary::MemPtr gathmethod;

  try {
    Arguments::const_iterator it = gathargs.lower_bound(str);
    Arguments::const_iterator to = gathargs.upper_bound(str);
    string pbc = "v";
    if (it == gathargs.upper_bound(str))
      throw Arguments::Exception("###### GATHER WARNING ######\n" + usage);
    else
      pbc = it->second;
    ++it;
    // define bool for checking
    bool list_available = false; // whether a list is available
    bool ref_available = false; // whether a ref is available
    bool uselist = false, useref = false; // conditioner

    if (it == gathargs.upper_bound(str)) {
      gathmethod = &Boundary::coggather;
      if(pbc != "v")
        std::cout << "###### GATHER WARNING ######\n"
	  	  << "# NO Gathering method specified ! \n"
		  << "# Thus : if the system is requested to be gathered, \n"
		  << "# the gathering will be done using the COG gather method\n";
    } else {
      std::string gather = it->second;
      if (writeout) cout << "# gather option : " << gather << endl;

      if (gather == "nog" || gather == "0") {
        gathmethod = &Boundary::nogather;
      } else if (gather == "1" || gather == "glist") {
        gathmethod = &Boundary::gatherlist;
        uselist = true;
      } else if (gather == "2" || gather == "gtime") {
        gathmethod = &Boundary::gathertime;
      } else if (gather == "3" || gather == "gref") {
        gathmethod = &Boundary::gatherref;
        useref = true;
      } else if (gather == "4" || gather == "gltime") {
        gathmethod = &Boundary::gatherltime;
        uselist = true;
      } else if (gather == "5" || gather == "grtime") {
        gathmethod = &Boundary::gatherrtime;
        useref = true;
      } else if (gather == "6" || gather == "gbond") {
        gathmethod = &Boundary::gatherbond;
      } else if (gather == "7" || gather == "cog") {
        gathmethod = &Boundary::coggather;
      } else if (gather == "8" || gather == "gfit") {
        gathmethod = &Boundary::gfitgather;
      } else {
        throw gromos::Exception("Gather", gather +
				" unknown. " + usage);
      }

      ++it;
      
      AtomSpecifier gathlist(sys);
      // first read all additional options (list, refg, molecules)
      if (it != gathargs.upper_bound(str)) {        
        while (it != to) {
          string gathopt = it->second;
          if (gathopt == "list") {
            it++;
            if (it == to) {
              throw gromos::Exception("Gathering ", "option " +gathopt +
				      ", but no list given!");
            }
            bool endlist=false;
            while (!endlist && it != to) {
              string gathopt = it->second.c_str();
              if (gathopt == "refg" || gathopt == "molecules") {
                endlist=true;
              } else {
                gathlist.addSpecifierStrict(gathopt);
                list_available=true;
                it++;
              }
            }
          } else if (gathopt == "refg") {
            // the reference frame itself is read in in the BoundaryParser
            it++;
            if (it == to) {
              throw gromos::Exception("Gathering ", "option " +gathopt +
				      ", but no reference structure given!");
            }
            ref_available=true;
            it++;
            
	  } else if (gathopt == "molecules") {
            // the specified molecules are read in in the BoundaryParser
	    it++;
	    if (it == to) {
	      throw gromos::Exception("Gathering ", "option " + gathopt +
				      ", but no molecules specified!");
	    }
	    it++;          
          } else {
            throw gromos::Exception("Unknown gather option ", gathopt +
				    " \nKnown gather options are : list, refg or molecules!");
          }
        }
      }
      
      if ((useref || gather == "gfit" || gather == "8") && list_available) {
        std::cout << "# Ignoring list option with gather method " << gather << "!" << std::endl;
      }
      if (uselist && ref_available) {
        std::cout << "# Ignoring refg option with gather method " << gather << "!" << std::endl;
      }
      if (gather == "5" || gather == "gtime") {
        if (ref_available && list_available) {
	  throw gromos::Exception("Gathering ",
				  "can not use both refg and list options at the same time, choose one!");
        }
        else if (ref_available) {
	  cout << "# NOTE: using gathermethod gtime with refg argument is equivalent to gather option grtime" << endl;
	  gathmethod = &Boundary::gatherrtime;
	  useref=true;
        }
        else if (list_available) {
	  cout << "# NOTE: using gathermethod gtime with list argument is equivalent to gather option gltime" << endl;
	  gathmethod = &Boundary::gatherltime;
	  uselist=true;
        }
      }
      if (uselist) {
        if (list_available) {
	  //the block of primlist
	  for (unsigned int j = 0; j < gathlist.size() / 2; ++j) {
	    int i = 2 * j;
	    sys.primlist[gathlist.mol(i)][0] = gathlist.atom(i);
	    sys.primlist[gathlist.mol(i)][1] = gathlist.mol(i + 1);
	    sys.primlist[gathlist.mol(i)][2] = gathlist.atom(i + 1);

	    refSys.primlist[gathlist.mol(i)][0] = gathlist.atom(i);
	    refSys.primlist[gathlist.mol(i)][1] = gathlist.mol(i + 1);
	    refSys.primlist[gathlist.mol(i)][2] = gathlist.atom(i + 1);

	    cout << "# updated prim : mol " << gathlist.mol(i) << " atom " << gathlist.atom(i)
		 << "\t# refe : mol " << sys.primlist[gathlist.mol(i)][1]
		 << " atom " << sys.primlist[gathlist.mol(i)][2] << endl;
	  }
        } else {        
	  cout << "# Gathering : You have requested to gather the system based on " << endl
	       << "# an atom list, while you didn't define such a list, therefore " << endl
	       << "# the gathering will be done according to the 1st atom of the previous molecule" << endl;
        }
      }
      if (useref && !ref_available) {
        throw gromos::Exception("Gathering ", gather +
				" No reference structure given!");
      }
    }
  } catch (Arguments::Exception &e) {
    gathmethod = &Boundary::gatherlist;
  }


  return gathmethod;
}

