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
 * @file boundary_checks.h
 * checks for boxes and boundary conditions
 */

#ifndef INCLUDED_BOUNDARY_CHECKS_H
#define	INCLUDED_BOUNDARY_CHECKS_H

namespace math {
  /**
   * checks whether the box is large enough for a given cutoff
   * @param box the box
   * @param b the boundary conditions
   * @param cutoff the cutoff
   */
  bool boundary_check_cutoff(math::Box const & box, math::boundary_enum const b,
          double cutoff);
}

#endif	/* _BOUNDARY_CHECKS_H */

