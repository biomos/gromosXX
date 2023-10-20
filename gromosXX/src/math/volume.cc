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
 * @file volume.cc
 */

#include "../stdheader.h"
#include "volume.h"

#undef MODULE
#undef SUBMODULE
#define MODULE math
#define SUBMODULE math

double math::volume(math::Box const & box, math::boundary_enum const b)
{
  switch (b){
    case math::vacuum:
      DEBUG(9, "vacuum: volume = 0");
      return 0;
    case math::rectangular:
      DEBUG(9, "rectangular: volume = "
	    <<  abs(box(0)) * abs(box(1)) * abs(box(2)));
      return abs(box(0)) * abs(box(1)) * abs(box(2));
    case math::truncoct:
    case math::triclinic:
      DEBUG(9, "triclinic: volume = "
	    << math::dot(math::cross(box(0), box(1)), box(2)));
      return math::dot(math::cross(box(0), box(1)), box(2));
    default:
      DEBUG(9, "volume error....");
  }
  return 0;
}
