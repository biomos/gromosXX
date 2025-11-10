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
 * @file InBFactorOccupancy.h
 * a file to read B-factor and occupancy data
 */

#ifndef INCLUDED_INBFACTOROCCUPANCY_H
#define	INCLUDED_INBFACTOROCCUPANCY_H

#include <string>
#include <vector>

#include "Ginstream.h"

namespace gio {
  /**
   * @struct BFactorOccupancyData
   * @ingroup gio
   * @author N.Schmid, F. Freitag
   *
   * A class to hold B-factor and occupancy data.
   */
  struct BFactorOccupancyData {
    /**
     * the B factor
     */
    double b_factor;
    /**
     * the occupancy value
     */
    double occupancy;
  };

  /**
   * @class InBFactorOccupancy
   * @ingroup gio
   * @author N. Schmid, F. Freitag
   * @brief reads a B-factor and Occupancy file
   * A class to read files containing B-factor and occupancy information
   *
   * Format of the B-factor and occupancy file:
   * The B-factors have to be given @f$\mathrm{nm}^2@f$.
   * @verbatim
TITLE
B-factors and occupancies for all atoms
END
BFACTOROCCUPANCY
# B-factor Occupancy
0.01  1.0
0.02  0.8
END
   @endverbatim
   */
  class InBFactorOccupancy {
  public:
    /**
     * constructor
     */
    InBFactorOccupancy() {}
    /**
     * construct from file name
     * @param file file name
     */
    InBFactorOccupancy(std::string file);
    /**
     * destructor
     */
    ~InBFactorOccupancy();
    /**
     * open a file
     * @param file file name
     */
    void open(std::string file);
    /**
     * close the file
     */
    void close();
    /**
     * get the mapping data
     * @return a map containing the IAC to element mapping
     */
    std::vector<BFactorOccupancyData> getData();
  protected:
    Ginstream file_stream;
  };
}

#endif	/* INCLUDED_INBFACTOROCCUPANCY_H */

