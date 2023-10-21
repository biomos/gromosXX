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
 * @file in_distanceres.h
 * read in a distrance restraints file.
 */
/**
 * @page disres distance restraints format
 * @date 28-10-2008
 *
 * A distance restraints specification file may contain the following blocks:
 * - @ref title
 * - @ref distanceresspec
 * - @ref pertdisresspec
 * - @ref mdisresspec
 * - @ref dfresspec
 * - @ref pertdfresspec
 */

#ifndef INCLUDED_IN_DISTANCERES_H
#define INCLUDED_IN_DISTANCERES_H

#include "../instream.h"

namespace io {

  /**
   * @class In_Distanceres
   * reads in a position restraints file
   */
  class In_Distanceres : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Distanceres() {}
    /**
     * Constructor.
     */
    In_Distanceres(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a distance restraints file.
     */
    void read(topology::Topology &topo,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);

    /**
     * read distance restraint specification block.
     */
    void read_DISTANCERESSPEC(topology::Topology &topo, simulation::Simulation &sim, std::ostream & os = std::cout);
    void read_PERTDISRESSPEC(topology::Topology &topo, simulation::Simulation &sim, std::ostream & os = std::cout);
    void read_DFRESSPEC(topology::Topology &topo, simulation::Simulation &sim, std::ostream & os = std::cout);
    void read_PERTDFRESSPEC(topology::Topology &topo, simulation::Simulation &sim, std::ostream & os = std::cout);
    void read_MDISRESSPEC(topology::Topology &topo, simulation::Simulation &sim, std::ostream & os = std::cout);

    /**
     * Maximum number of atoms that can be specified to define a virtual atom
     */
    static const unsigned int MAX_ATOMS = 4;
  };

} // io

#endif
