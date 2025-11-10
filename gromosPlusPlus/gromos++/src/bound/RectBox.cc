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

// bound_RectBox.cc

#include <cmath>

#include "RectBox.h"
#include "../gmath/Vec.h"
#include "../gcore/Box.h"

using bound::RectBox;
using gmath::Vec;
using gcore::Box;

Vec RectBox::nearestImage(const Vec &r1, const Vec &r2, const Box &box)const {
  Vec diff = r2 - r1;
  Vec a;
  const double kabs = box.K().abs();
  a[0] = diff[0] - kabs * rint(diff[0] / kabs);
  const double labs = box.L().abs();
  a[1] = diff[1] - labs * rint(diff[1] / labs);
  const double mabs = box.M().abs();
  a[2] = diff[2] - mabs * rint(diff[2] / mabs);

  Vec rec = r1 + a;
  return rec;
}


