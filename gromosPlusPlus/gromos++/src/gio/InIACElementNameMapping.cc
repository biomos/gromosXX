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
 * @file InIACElementNameMapping.cc
 * implemenation of the element mapping reader
 */

#include "InIACElementNameMapping.h"

#include <cstddef>
#include <map>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

#include "../gromos/Exception.h"

gio::InIACElementNameMapping::InIACElementNameMapping(std::string file) {
  open(file);
}

gio::InIACElementNameMapping::~InIACElementNameMapping() {
  close();
}

void gio::InIACElementNameMapping::open(std::string file) {
  file_stream.open(file);
}

void gio::InIACElementNameMapping::close() {
  file_stream.close();
}

std::map<int,std::string> gio::InIACElementNameMapping::getData() {
  std::vector<std::string> buffer;
  std::map<int, std::string> gactoele;

  file_stream.getblock(buffer);
  if (buffer[0] != "ELEMENTMAPPING")
    throw gromos::Exception("getElementMapping",
          "library file does not contain a ELEMENTMAPPING block!");
  if (buffer[buffer.size() - 1].find("END") != 0)
    throw gromos::Exception("InIACElementNameMapping", "No END in ELEMENTMAPPING"
          " block.");
  // read in the lib
  for (size_t i = 1; i < buffer.size() - 1; i++) {
    int gac;
    std::string element;
    std::istringstream ss(buffer[i]);
    ss >> gac >> element;
    if (ss.fail())
      throw gromos::Exception("InIACElementNameMapping", "bad line in ELEMENTMAPPING"
            " block.");
    gactoele[gac - 1] = element;
  }
  return gactoele;
}



