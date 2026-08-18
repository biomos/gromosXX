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

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

#include "tre_parser.h"

namespace {
  // Trims trailing '\r' (harmless either way) and surrounding
  // whitespace for the block-marker string comparisons below.
  std::string trimmed(const std::string & s) {
    size_t b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    size_t e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
  }
}

bool check_ubiquitin::parse_tre(const std::string & path,
                                 std::vector<check_ubiquitin::TreStep> & out,
                                 std::ostream & err) {
  std::ifstream f(path);
  if (!f) {
    err << "tre_parser: could not open " << path << std::endl;
    return false;
  }
  out.clear();

  std::string line;
  bool have_step = false;
  unsigned current_step = 0;
  double current_time = 0.0;

  while (std::getline(f, line)) {
    const std::string t = trimmed(line);
    if (t == "TIMESTEP") {
      if (!std::getline(f, line)) {
        err << "tre_parser: " << path << ": truncated TIMESTEP block" << std::endl;
        return false;
      }
      std::istringstream iss(line);
      if (!(iss >> current_step >> current_time)) {
        err << "tre_parser: " << path << ": malformed TIMESTEP line: " << line
            << std::endl;
        return false;
      }
      have_step = true;
      if (!std::getline(f, line) || trimmed(line) != "END") {
        err << "tre_parser: " << path << ": TIMESTEP block missing END" << std::endl;
        return false;
      }
    } else if (t == "ENERGY03") {
      if (!have_step) {
        err << "tre_parser: " << path << ": ENERGY03 block with no preceding "
               "TIMESTEP" << std::endl;
        return false;
      }
      if (!std::getline(f, line) || trimmed(line) != "# totals") {
        err << "tre_parser: " << path << ": ENERGY03 block missing \"# totals\""
            << std::endl;
        return false;
      }
      TreStep step;
      step.step = current_step;
      step.time = current_time;
      for (int i = 0; i < TOTALS_FIELD_COUNT; ++i) {
        if (!std::getline(f, line)) {
          err << "tre_parser: " << path << ": truncated \"# totals\" block ("
              << i << "/" << TOTALS_FIELD_COUNT << " fields read)" << std::endl;
          return false;
        }
        std::istringstream iss(line);
        if (!(iss >> step.totals[i])) {
          err << "tre_parser: " << path << ": non-numeric totals field: " << line
              << std::endl;
          return false;
        }
      }
      // Per-energy-group-pair matrices etc. follow -- not needed here,
      // just consume lines until this block's own END.
      while (std::getline(f, line) && trimmed(line) != "END") {}
      out.push_back(step);
    } else if (t == "VOLUMEPRESSURE03") {
      if (out.empty()) {
        err << "tre_parser: " << path << ": VOLUMEPRESSURE03 block with no "
               "preceding ENERGY03" << std::endl;
        return false;
      }
      while (std::getline(f, line) && trimmed(line) != "# temperature") {}
      unsigned num_baths = 0;
      if (!std::getline(f, line)) {
        err << "tre_parser: " << path << ": truncated VOLUMEPRESSURE03 block"
            << std::endl;
        return false;
      }
      {
        std::istringstream iss(line);
        if (!(iss >> num_baths)) {
          err << "tre_parser: " << path << ": malformed bath count: " << line
              << std::endl;
          return false;
        }
      }
      std::vector<TreBath> baths(num_baths);
      for (unsigned b = 0; b < num_baths; ++b) {
        if (!std::getline(f, line)) {
          err << "tre_parser: " << path << ": truncated bath temperature line"
              << std::endl;
          return false;
        }
        std::istringstream iss(line);
        if (!(iss >> baths[b].total >> baths[b].com >> baths[b].ir >> baths[b].scale)) {
          err << "tre_parser: " << path << ": malformed bath temperature line: "
              << line << std::endl;
          return false;
        }
      }
      out.back().baths = baths;
      while (std::getline(f, line) && trimmed(line) != "END") {}
      have_step = false;
    }
  }

  if (out.empty()) {
    err << "tre_parser: " << path << ": no ENERGY03 blocks found" << std::endl;
    return false;
  }
  return true;
}

bool check_ubiquitin::scan_trc_positions(const std::string & path, std::ostream & err) {
  std::ifstream f(path);
  if (!f) {
    err << "tre_parser: could not open " << path << std::endl;
    return false;
  }

  constexpr double kMaxAbsCoord = 1.0e4; // nm -- generous, real coords here are ~O(1-10)

  std::string line;
  long block_count = 0;
  long line_no = 0;
  bool in_block = false;
  while (std::getline(f, line)) {
    ++line_no;
    const std::string t = trimmed(line);
    if (t == "POSITIONRED") {
      in_block = true;
      ++block_count;
      continue;
    }
    if (!in_block) continue;
    if (t == "END") {
      in_block = false;
      continue;
    }
    if (t.empty() || t[0] == '#') continue;

    std::istringstream iss(line);
    double x, y, z;
    if (!(iss >> x >> y >> z)) {
      err << "tre_parser: " << path << ":" << line_no
          << ": malformed POSITIONRED line: " << line << std::endl;
      return false;
    }
    for (double v : {x, y, z}) {
      if (std::isnan(v) || std::isinf(v) || std::abs(v) > kMaxAbsCoord) {
        err << "tre_parser: " << path << ":" << line_no
            << ": position out of sane range: " << line << std::endl;
        return false;
      }
    }
  }

  if (block_count == 0) {
    err << "tre_parser: " << path << ": no POSITIONRED blocks found" << std::endl;
    return false;
  }
  return true;
}
