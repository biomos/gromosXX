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

#ifndef INCLUDED_UTILS_GCH
#define INCLUDED_UTILS_GCH
#include <vector>
#include <string>

#include "../gromos/Exception.h"
#include "../gcore/System.h"

namespace gcore{
  class System;
  class GromosForcefield;
  class Bond;
  class Angle;
}

namespace gmath{
  class Vec;
}

namespace utils
{
  /**
  * fill the hydrogen and non-hydrogen neighbours of atom a of molecule m
  * into vectors h and nh
  */
  void get_h_nh_neighbours(gcore::System &sys, gcore::GromosForceField &gff, int m, int a, std::vector<int> &h, std::vector<int> &nh);
   
  /**
  * determine the geometry of an atom based on the number of hydrogen
  * and non-hydrogen neighbours
  */ 
  int get_geometry(int numH, int numNH);

  /**
  * recalculate optimal positions for hydrogens based on the 
  * @param m molecule number
  * @param a atom number
  * @param h vector of hydrogen neighbours of a
  * @param nh vector of non-hydrogen neighbours of a
  * @param geom geometry type, see @ref gch 
  * @param eps tolerated deviation from optimal bond length
  */
  int generate_hcoordinates(gcore::System &sys, gcore::GromosForceField &gff, int m, int a, std::vector<int> &h, std::vector<int> &nh, int geom, double eps);  
  int generate_hcoordinates(gcore::System &sys, gcore::GromosForceField &gff, int m, int a, double eps);

}
#endif
