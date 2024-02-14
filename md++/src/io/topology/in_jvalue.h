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
 * @file in_jvalue.h
 * read in a jvalue restraining specification file
 */
/**
 * @page jvalue J-value restraints specification format
 * @date 28-10-2008
 *
 * A J-value restraining specifcation file may contain the following
 * blocks:
 * - @ref title
 * - @ref jvalresspec
 */

#ifndef INCLUDED_IN_JVALUE_H
#define INCLUDED_IN_JVALUE_H

#include "../instream.h"

namespace io {

  /**
   * @class In_Jvalue
   * Reads in a @f$^3J@f$-value restraining specification file.
   */
  class In_Jvalue : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Jvalue() {}
    /**
     * Constructor.
     */
    In_Jvalue(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a @f$^3J@f$-value restraining specification file.
     */
    void read(topology::Topology &topo,
	      configuration::Configuration & conf,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);

  };

} // io

#endif
