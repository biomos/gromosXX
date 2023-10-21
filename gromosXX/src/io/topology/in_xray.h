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
 * @file in_xray.h
 * read in a x-ray restraining file.
 */
/**
 * @page xrayresfile X-ray restraints format
 * @date 12-07-2011
 *
 * A x-ray restraints specification file may contain the following
 * blocks:
 * - @ref title
 * - @ref xrayresspec
 * - @ref xrayresspec "XRAYRFREESPEC block"
 * - @ref elementmap
 * - @ref elementmap "XRAYSOLVELEMENTSPEC block"
 * - @ref xraybfoccspec
 * - @ref xraybfoccspec "XRAYSOLVBFOCCSPEC block"
 * - @ref xrayrespara
 * - @ref xrayumbrellaweight
 * - @ref xraysymresspec
 * - @ref xraybfopt
 * - @ref xraysfcalc
 * - @ref xrayreplicaexchange
 * - @ref xrayoverallbfactor
 */

#ifndef INCLUDED_IN_XRAYRES_H
#define INCLUDED_IN_XRAYRES_H

#include "../instream.h"

namespace io {

  /**
   * @class In_Xrayresspec
   * reads in a xray restraining specification file
   */
  class In_Xrayresspec : public GInStream {

  public:
    /**
     * Default constructor.
     */
    In_Xrayresspec() {}
    /**
     * Constructor.
     */
    In_Xrayresspec(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a xray restraining specification file.
     */
    void read(topology::Topology &topo,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);

  };
  
} // io

#endif


