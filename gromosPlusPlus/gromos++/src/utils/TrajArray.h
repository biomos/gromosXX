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

// TrajArray.h
#ifndef TRAJARRAY_H
#define TRAJARRAY_H

#include "../gcore/Box.h"
#include "../gcore/Molecule.h"
#include "../gcore/System.h"

namespace utils {
  /**
   * Class TrajArray
   * The TrajArray contains a complete trajectory
   *
   * If the coordinates of multiple frames is needed, the TrajArray can
   * be used to access these simultaneously without needing to store
   * the topological information of your system hundreds of times
   *
   * @class TrajArray
   * @author T. Hansson and V. Kraeutler
   * @ingroup utils
   */
  class TrajArray {
  public:
    /**
     * Constructor
     */
    TrajArray(const gcore::System &sys);

    /**
     * Destructor
     **/
    ~TrajArray();

    /**
     * store a frame in the array
     * @param sys the system to store
     * @param frameIndex index of the frame
     */
    void store(const gcore::System &sys,
        const unsigned int frameIndex);

    /**
     * extract a system from the array
     * @param sys the resulting system
     * @param frameIndex the index of the frame you want to extract
     */
    void extract(gcore::System &sys,
        const unsigned int frameIndex) const;

    /**
     * the number of of coordinates stored per frame
     */
    inline unsigned int numAtoms();

  protected:
    std::vector<double *> trajectoryData;
    // number of coordinates per frame
    unsigned int nAtoms;

  };
}
#endif                                                                 
