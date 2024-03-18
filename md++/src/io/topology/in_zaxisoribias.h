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
 * @file in_zaxisoribias.h
 * read in a z-axis orientation bias file.
 */
/**
 * @page zaxisori z-axis orientation bias format
 * @date 18-06-2019
 *
 * A z-axis orientation bias specification file may contain the following blocks:
 * - @ref title
 * - @ref zaxisoribiasspec
 */

#ifndef INCLUDED_IN_ZAXISORIBIAS_H
#define INCLUDED_IN_ZAXISORIBIAS_H

#include "../instream.h"

namespace io {

  /**
   * @class In_Zaxisoribias
   * reads in a z-axis orientation bias file
   */
  class In_Zaxisoribias : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Zaxisoribias() {}
    /**
     * Constructor.
     */
    In_Zaxisoribias(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a z-axis orientation bias file.
     */
    void read(topology::Topology &topo,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);

    /**
     * read z-axis orientation bias specification block.
     */
    void read_ZAXISORIBIASSPEC(topology::Topology &topo, simulation::Simulation &sim, std::ostream & os = std::cout);

    /**
     * Maximum number of atoms that can be specified to define a virtual atom
     */
    static const unsigned int MAX_ATOMS = 4;
  };

} // io

#endif
