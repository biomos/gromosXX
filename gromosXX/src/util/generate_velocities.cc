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
 * @file generate_velocities.cc
 * generate velocities from Maxwell distribution.
 */

#include <stdheader.h>
#include <math/random.h>
#include <simulation/multibath.h>
#include <simulation/parameter.h>

#include "generate_velocities.h"

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE util

void util::generate_velocities(const simulation::Parameter &param,
                               double const temp, math::SArray const &mass,
			       math::VArray & vel, math::VArray & old_vel,
			       unsigned int const seed,
			       std::ostream & os, bool quiet)
{

  if (!quiet) {
    os << "\n\tgenerating initial velocities\n"
       << "\t\ttemperature         = " << temp << "\n"
       << "\t\trandom number seed  = " << seed << "\n";
  }
  
  std::ostringstream stringseed; stringseed << seed;
  math::RandomGenerator* rng = math::RandomGenerator::create(param, stringseed.str());
  if (!quiet)
    os << "\t\trandom number generator: " << rng->description() << "\n";
  
  for(unsigned int i=0; i<vel.size(); ++i){
    const double sd = sqrt(math::k_Boltzmann * temp / mass(i));
    rng->stddev(sd);
    for(int d=0; d<3; ++d){
      // old_vel(i)(d) = sd * n.random();
      old_vel(i)(d) = rng->get_gauss();
      vel(i)(d) = old_vel(i)(d);
    }
  }

  delete rng;
}
