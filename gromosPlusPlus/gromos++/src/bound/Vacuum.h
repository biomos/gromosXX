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

// bound_Vacuum.h

#ifndef INCLUDED_BOUND_VACUUM
#define INCLUDED_BOUND_VACUUM


#include "Boundary.h"
#include "../gmath/Vec.h"

namespace bound {

  /**
   * Class Vacuum
   * defines the periodic boundary conditions for a vacuum. Which means 
   * that there are no periodic boundary conditions
   *
   * @class Vacuum
   * @author R. Buergi
   * @ingroup bound
   */
  class Vacuum : public Boundary {
  public:
    Vacuum(gcore::System *sys) : Boundary(sys) {
      setType('v');
    }

    virtual ~Vacuum() {
    }

    virtual gmath::Vec nearestImage(const gmath::Vec &r1,
            const gmath::Vec &r2,
            const gcore::Box &box) const {
      return r2;
    }

    // overwrite gathering methods as they do not make sense for vacuum
    virtual void nogather() {
    }

    virtual void gathergr() {
    }

    virtual void gather() {
    }

    virtual void coggather() {
    }

    virtual void seqgather() {
    }

    virtual void crsgather() {
    }

    virtual void gengather() {
    }

    virtual void gatherlist() {
    }

    virtual void gathertime() {
    }

    virtual void gatherref() {
    }

    virtual void gatherltime() {
    }

    virtual void gatherrtime() {
    }

    virtual void gatherbond() {
    }
  };

}

#endif
