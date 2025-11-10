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

// args_BoundaryParser.cc
#include "BoundaryParser.h"

#include <string>
#include <vector>

#include "Arguments.h"
#include "../bound/Vacuum.h"
#include "../bound/TruncOct.h"
#include "../bound/RectBox.h"
#include "../bound/Triclinic.h"
#include "../gcore/System.h"
#include "../gcore/Box.h"
#include "../utils/parse.h"
#include "../gromos/Exception.h"

using namespace args;
using namespace bound;

bound::Boundary *BoundaryParser::boundary(gcore::System &sys,
        const Arguments &args,
        const std::string &str) {
  Boundary *pbc;

  try {

    // get the boundary from a read box, only if no @pbc flag is set in the args
    if (sys.hasBox && sys.box().boxformat() != gcore::Box::box96 && args.count("pbc") < 0) {
      switch (sys.box().ntb()) {
        case gcore::Box::vacuum:
          pbc = new Vacuum(&sys);
          break;
        case gcore::Box::rectangular:
          pbc = new RectBox(&sys);
          break;
        case gcore::Box::triclinic:
          pbc = new Triclinic(&sys);
          break;
        case gcore::Box::truncoct:
          pbc = new TruncOct(&sys);
          break;
        default:
          throw gromos::Exception("Boundary", "Could not parse the boundary from"
                  " the box.");
      }
    } else { // parse the boundary condition with help of the command line argument @pbc
      Arguments::const_iterator it = args.lower_bound(str);
      if (it == args.upper_bound(str))
        throw Arguments::Exception("");

      switch (it->second[0]) {
        case 'r':
          pbc = new RectBox(&sys);
          break;
        case 't':
          pbc = new TruncOct(&sys);
          break;
        case 'c':
          pbc = new Triclinic(&sys);
          break;
        case 'v':
          pbc = new Vacuum(&sys);
          break;
        default:
          throw gromos::Exception("Boundary", args[str] +
                  " unknown. Known boundaries are r (rectangular), "
                  "t (truncated octahedron), c (triclinic) and v (vacuum)");
      }

      ++it;
    }
  } catch (Arguments::Exception &e) {
    pbc = new Vacuum(&sys);
  }

  // this is for additional gathering flags
  Arguments::const_iterator it = args.lower_bound(str),
          to = args.upper_bound(str);
  std::string gmethod = "";
  if (it != to) {
    ++it; // that was the boxshape
    if (it != to) {
      gmethod = it->second;
      it++;
    }
  }
  bool mol=false;
  bool ref=false;
  // reference frame is the first frame by default, overwritten below by the refg argument if present
  std::string refframe;
  for (; it != to; ++it) {
    if (it->second == "refg") { 
      ref=true;
      ++it;
      if (it == to) {
        throw Arguments::Exception("You have to specify a filename when using refg gathering method");
      }
      pbc->setReferenceFrame(it->second);
    }
    else if (it->second == "molecules" && gmethod == "gfit") {
      mol=true;
      ++it;
        std::vector<int> molecules;
        utils::parse_range(it->second.c_str(), molecules);
        for (unsigned int i=0; i<molecules.size(); i++) {
          pbc->addRefMol(molecules[i]);
          //std::cout << "# Gathering molecule " << molecules[i] << " to reference" << std::endl;
        }
    }
  }
  if (!ref && (gmethod == "gfit" || gmethod == "2" ||gmethod == "gtime" || gmethod == "2" )) {
    pbc->setReferenceFrame(args.lower_bound("traj")->second);
  }
  
  // if no molecules are given to gather to reference
  // do it with all solute molecules
  if (!mol) {
    //std::cout << "# no molecules selected, by default gathering all solute molecules to reference in first step" << std::endl;
    for (int i = 1; i-1 < sys.numMolecules(); ++i) {
      pbc->addRefMol(i);
    }
  }
  return pbc;
}

bound::Boundary *BoundaryParser::boundary(gcore::System &sys) {

  Boundary *pbc;
  // get the boundary from a read box, only if no @pbc flag is set in the args
  if (sys.hasBox && sys.box().boxformat() != gcore::Box::box96) {
    switch (sys.box().ntb()) {
      case gcore::Box::vacuum:
        pbc = new Vacuum(&sys);
        break;
      case gcore::Box::rectangular:
        pbc = new RectBox(&sys);
        break;
      case gcore::Box::triclinic:
        pbc = new Triclinic(&sys);
        break;
      case gcore::Box::truncoct:
        pbc = new TruncOct(&sys);
        break;
      default:
        throw gromos::Exception("Boundary", "Could not parse the boundary from"
                " the box.");
    }
  } else {
    throw gromos::Exception("Boundary", "Could not parse the boundary from the box.");
  }

  return pbc;
}
