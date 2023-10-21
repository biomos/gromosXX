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
 * @file in_posres.h
 * read in a position restraining file.
 */
/**
 * @page posres position restraints format
 * @date 28-10-2008
 *
 * A position restraints/constraints specifcation file may contain the following
 * blocks:
 * - @ref title
 * - @ref posresspec
 * - @ref refposition
 * - @ref bfactor
 */

#ifndef INCLUDED_IN_POSRES_H
#define INCLUDED_IN_POSRES_H

#include "../instream.h"

namespace io {

  /**
   * @class In_Posresspec
   * reads in a position restraining specification file
   */
  class In_Posresspec : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Posresspec() {}
    /**
     * Constructor.
     */
    In_Posresspec(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a position restraining specification file.
     */
    void read(topology::Topology &topo,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);

  };
  
  /**
   * @class In_Refpos
   * reads in a position restraining file
   */
  class In_Refpos : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Refpos() {}
    /**
     * Constructor.
     */
    In_Refpos(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a position restraining file.
     */
    void read(topology::Topology &topo,
            configuration::Configuration & conf,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);

  };

} // io

#endif
