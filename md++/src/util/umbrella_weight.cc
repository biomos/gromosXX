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
 * @file umbrella_weight.cc
 * implementation of umbrella weighting
 */
#include "../stdheader.h"

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"
#include "../util/debug.h"

#include "umbrella_weight.h"

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE leus

std::ostream & operator<<(std::ostream & os, const util::Umbrella_Weight & w) {
  w.write(os);
  return os;
}
std::istream & operator>>(std::istream & is, util::Umbrella_Weight & w) {
  w.read(is);
  return is;
}

void util::Number_Of_Visits_Umbrella_Weight::write(std::ostream & os) const {
  os << std::setw(10) << weight;
}

util::Umbrella_Weight_Factory::~Umbrella_Weight_Factory() {
  // delete all the created instaces
  for(std::vector<util::Umbrella_Weight*>::iterator it = instances.begin(),
          to = instances.end(); it != to; ++it)
    delete *it;
}

