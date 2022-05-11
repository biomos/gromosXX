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
 * @file in_rdc.h
 * read in a rdc restraining specification file.
 */
/**
 * @page rdc RDC restraints specification format
 * @date 09-05-2011
 *
 * A RDC restraints specification file has to contain the following
 * blocks:
 * - @ref title
 * - @ref conversion
 * - @ref magfieldc
 * - @ref alignt
 * - @ref rdcresspec
 */

#ifndef INCLUDED_IN_RDC_H
#define INCLUDED_IN_RDC_H

#include "../instream.h"

namespace io {

  /**
   * @class In_RDC
   * reads in a RDC restraint specification file.
   */
  class In_RDC : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_RDC() {}
    /**
     * Constructor.
     */
    In_RDC(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a RDC restraining file.
     */
    void read(topology::Topology &topo,
	      configuration::Configuration & conf,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);
    /**
     * Maximum number of atoms that can be specified to define a virtual atom
     */
    static const unsigned int MAX_ATOMS = 4;

  };

} // io

#endif

