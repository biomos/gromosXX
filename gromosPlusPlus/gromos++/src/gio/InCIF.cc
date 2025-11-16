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
 * @file InCIF.cc
 * implemenation of the CIF reader
 */

#include "InCIF.h"

#include <sstream>
#include <iostream>
#include <string>
#include <vector>

#include "StringTokenizer.h"
#include "../gromos/Exception.h"

gio::InCIF::InCIF(std::string file) {
  open(file);
}

gio::InCIF::~InCIF() {
  if (isOpen()) close();
}

void gio::InCIF::open(std::string file) {
  file_stream.open(file.c_str());
  if (!file_stream.is_open())
    throw gromos::Exception("InCIF", "Cannot open CIF file " + file + ".");
}

void gio::InCIF::close() {
  if (isOpen())
    file_stream.close();
}

bool gio::InCIF::isOpen() {
  return file_stream.is_open();
}

std::vector<gio::CIFData> gio::InCIF::getData() {
  std::vector<gio::CIFData> cifdata;
  gio::CIFData tcif;
  std::string line;
  int refln_nr = 0, hnr = 0, knr = 0, lnr = 0, Fnr = 0, Fsnr = 0;
  while (true) {
    if (file_stream.eof()) {
      break;
    }
    std::getline(file_stream, line);
    if (line.substr(0, 6) == "_refln") {
      refln_nr++;
      // check index of desired data
      if (line == "_refln.index_h") {
        hnr = refln_nr;
      } else if (line == "_refln.index_k") {
        knr = refln_nr;
      } else if (line == "_refln.index_l") {
        lnr = refln_nr;
      } else if (line == "_refln.F_meas_au") {
        Fnr = refln_nr;
      } else if (line == "_refln.F_meas_sigma_au") {
        Fsnr = refln_nr;
      }
    } else {
      // Check if before or after _refln-Block
      if (refln_nr == 0) {
        // do nothing / del lines
      } else {
        // check if data or not
        gio::StringTokenizer tok(line);
        std::vector<std::string> split = tok.tokenize();
        if (split.size() < 4) {
          // no valid line...skip
        } else {
          std::istringstream is1(split[hnr - 1]); is1 >> tcif.h;
          std::istringstream is2(split[knr - 1]); is2 >> tcif.k;
          std::istringstream is3(split[lnr - 1]); is3 >> tcif.l;
          std::istringstream is4(split[Fnr - 1]); is4 >> tcif.f_obs;
          if (Fsnr != 0) {
            std::istringstream is(split[Fsnr - 1]); is >> tcif.stddev_f_obs;
          } else {
            tcif.stddev_f_obs = 0.0;
          }
          if (is1.fail() || is2.fail() || is3.fail() || is4.fail()) 
            std::cerr << "Skipping bad line." << std::endl;
          else
            cifdata.push_back(tcif);
        }
      }
    }
  }
  return cifdata;
}


