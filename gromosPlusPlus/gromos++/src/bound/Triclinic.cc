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

// bound_Triclinic.cc

#include "Triclinic.h"

#include <cmath>

#include "../gmath/Vec.h"
#include "../gcore/Box.h"

using bound::Triclinic;
using gmath::Vec;
using gcore::Box;

Vec Triclinic::nearestImage(const Vec &r1, const Vec &r2, const Box &box)const{
  Vec P = r2 - r1;
  int k,l,m;
  k = int(rint(box.cross_K_L_M()[0].dot(P)));
  l = int(rint(box.cross_K_L_M()[1].dot(P)));
  m = int(rint(box.cross_K_L_M()[2].dot(P)));
  
  P += box.K() * k + box.L() * l + box.M() * m;
  return r1 + P;
}
