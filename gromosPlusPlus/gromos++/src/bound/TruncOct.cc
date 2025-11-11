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

// bound_TruncOct.cc

#include "TruncOct.h"

#include <cmath>

#include "../gmath/Vec.h"
#include "../gcore/System.h"
#include "../gcore/Solvent.h"
#include "../gcore/Molecule.h"
#include "../gcore/Box.h"

using bound::TruncOct;
using gmath::Vec;
using gcore::Box;

Vec TruncOct::nearestImage(const Vec &r1, const Vec &r2, const Box &box)const{
  Vec diff=r2-r1;
  Vec a;

  const double kabs = box.K().abs();
  a[0] = diff[0] - kabs * rint(diff[0]/kabs);
  a[1] = diff[1] - kabs * rint(diff[1]/kabs);
  a[2] = diff[2] - kabs * rint(diff[2]/kabs);

  if ( (0.75*kabs - fabs(a[0]) - fabs(a[1]) - fabs(a[2])) < 0.0) {
    const double half_kabs = 0.5 * kabs;
    a[0] = a[0] - a[0]/fabs(a[0])*half_kabs;
    a[1] = a[1] - a[1]/fabs(a[1])*half_kabs;
    a[2] = a[2] - a[2]/fabs(a[2])*half_kabs;

  }

  return r1 + a;
}
