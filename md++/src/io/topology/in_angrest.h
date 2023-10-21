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
 * @file in_angrest.h
 * read in a angle restraints file.
 */
/**
 * @page angrest angle restraints format
 * @date 24-05-2019
 *
 * A angle restraints/constraints specification file may contain the following
 * blocks:
 * - @ref title
 * - @ref angresspec
 * - @ref pertangresspec
 */
#ifndef INCLUDED_IN_ANGREST_H
#define INCLUDED_IN_ANGREST_H

#include "../instream.h"

namespace io {

  /**
   * @class In_Angrest
   * reads in a angle restraints file
   */
  class In_Angrest : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Angrest() {}
    /**
     * Constructor.
     */
    In_Angrest(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a position restraints file.
     */
    void read(topology::Topology &topo,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);

    void read_ANGRESSPEC(topology::Topology &topo,
        									simulation::Simulation &sim,
        									std::ostream & os);
    void read_PERTANGRESSPEC(topology::Topology &topo,
									simulation::Simulation &sim,
									std::ostream & os);
  };
} // io

#endif
