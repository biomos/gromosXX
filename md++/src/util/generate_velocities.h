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
 * @file generate_velocities.h
 * generate initial velocities from Maxwell distribution.
 */

#ifndef INCLUDED_GENERATE_VELOCITIES_H
#define INCLUDED_GENERATE_VELOCITIES_H

namespace simulation {
  class Parameter;
}

namespace util
{
  /**
   * generate initial velocities.
   */
  void generate_velocities(const simulation::Parameter &param,
                           double const temp, math::SArray const &mass,
			   math::VArray &vel, math::VArray &old_vel,
			   unsigned int const seed,
			   std::ostream & os = std::cout, 
                           bool quiet = false);
  
}

#endif