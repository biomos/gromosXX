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
 * @file InBFactorOccupancy.cc
 * implemenation of the b-factor and occupancy reader
 */

#include "InBFactorOccupancy.h"

#include <cstddef>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

#include "../gromos/Exception.h"

gio::InBFactorOccupancy::InBFactorOccupancy(std::string file) {
  open(file);
}

gio::InBFactorOccupancy::~InBFactorOccupancy() {
  close();
}

void gio::InBFactorOccupancy::open(std::string file) {
  file_stream.open(file);
}

void gio::InBFactorOccupancy::close() {
  file_stream.close();
}

std::vector<gio::BFactorOccupancyData> gio::InBFactorOccupancy::getData() {
  std::vector<std::string> buffer;
  std::vector<gio::BFactorOccupancyData> bfoccu;
  gio::BFactorOccupancyData sdata;
  // Get Block
  file_stream.getblock(buffer);
  if (buffer[0] != "BFACTOROCCUPANCY")
    throw gromos::Exception("InBFactorOccupancy",
          "library file does not contain a BFACTOROCCUPANCY block!");
  if (buffer[buffer.size() - 1].find("END") != 0)
    throw gromos::Exception("InBFactorOccupancy", "No END in BFACTOROCCUPANCY"
          " block.");
  // read in the lib
  for (size_t i = 1; i < buffer.size() - 1; i++) {
    // Data-Vars
    std::istringstream ss(buffer[i]);
    ss >> sdata.b_factor >> sdata.occupancy;
    if (ss.fail())
      throw gromos::Exception("InBFactorOccupancy", "bad line in BFACTOROCCUPANCY"
            " block.");
    // Push back to vetor
    bfoccu.push_back(sdata);
  }
  return bfoccu;
}




