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
 * @file in_symrest.h
 * read in a symmetry restraints specification file
 */
/**
 * @page symrest symmetry restraints specification format
 * @date 17-08-2011
 *
 * A symmetry restraints specification file may contain the following
 * blocks:
 * - @ref title
 * - @ref transform
 * - @ref symresspec
 */

#ifndef INCLUDED_IN_SYMREST_H
#define INCLUDED_IN_SYMREST_H

#include "../instream.h"

namespace io {

  /**
   * @class In_Symrest
   * reads in a symmetry restraints specification file
   */
  class In_Symrest : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Symrest() {}
    /**
     * Constructor.
     */
    In_Symrest(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a symmetry restraints specification file.
     */
    void read(topology::Topology &topo,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);

  };
} // io

#endif
