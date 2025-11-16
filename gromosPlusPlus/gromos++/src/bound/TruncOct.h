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

// bound_TruncOct.h

#ifndef INCLUDED_BOUND_TRUNCOCT
#define INCLUDED_BOUND_TRUNCOCT

#include "Boundary.h"

#include "../gromos/Exception.h"

namespace bound{
  /**
   * Class TruncOct
   * Class that defines periodic boundary conditions for a truncated 
   * octahedral box
   *
   * @class TruncOct
   * @author R. Buergi
   * @ingroup bound
   */
  class TruncOct: public Boundary {
  public:
    // Constructor
    TruncOct(gcore::System *sys): Boundary(sys){setType('t');}
    virtual ~TruncOct(){}
    virtual gmath::Vec nearestImage(const gmath::Vec &r1,
			    const  gmath::Vec &r2, 
			    const gcore::Box &box) const;
  };
    
}

#endif
