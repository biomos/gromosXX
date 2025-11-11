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

// args_GatherParser.h

#ifndef INCLUDED_ARGS_GATHERPARSER
#define INCLUDED_ARGS_GATHERPARSER

#include "../bound/Boundary.h"

#include <string>


namespace gcore{
  class System;
}

namespace bound{
  class Boundary;
}

namespace args{
  class Arguments;


/**
 * Class GatherParser
 * Purpose: Parse gathering methods from args.
 * 
 *
 * Description:
 * This class is used to parse gathering methods from args.
 * By default (when no arguments are given) it will use the gather-gromos 
 * gathering method. By default it will use args["pbc"].
 *
 * 
 * @class GatherParser
 * @version $Date: Mon Jul 15 14:17:41 MEST 2002
 * @author  M.A. Kastenholz
 * @ingroup args
 * @sa args::Arguments
 * @sa args::BoundaryParser
 */
  class GatherParser{
    
    // not implemented
    GatherParser(const GatherParser&);
    GatherParser &operator=(const GatherParser&);

  public:
/**
 * GatherParser constructor.
 * Details.
 */
    GatherParser(){}
/**
 * GatherParser destructor.
 * Details.
 */
    ~GatherParser(){}
/** 
 * Constructs the class and returns a member pointer to a gathering method.
 * Method parse parses the input from args.
 * @param gathargs Arguments from the input line.
 * @param str name of the argument string (default "pbc")
 * @return bound::Boundary::MemPtr Member pointer to the gathering method.
 * Details.
 */
    static bound::Boundary::MemPtr parse(gcore::System &sys,gcore::System &refSys,const Arguments &gathargs, const std::string &str = "pbc", const bool writeout=0);

  };

}

#endif
