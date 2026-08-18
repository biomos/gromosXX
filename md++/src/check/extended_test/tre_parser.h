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
 * @file tre_parser.h
 * Minimal reader for the .tre energy trajectory format
 * (io::Out_Configuration::_print_energyred_helper,
 * src/io/configuration/out_configuration.cc) -- exists specifically so
 * tests can check the *actual bytes program/md writes to disk*, not
 * just the in-memory configuration::Energy values a test harness reads
 * directly. Those are not the same code path: a bug in
 * Out_Configuration's own formatting/precision, or in which
 * configuration object it reads from, would be invisible to a test
 * that only inspects memory (this is exactly why ubiquitin_cpu.t.cc/
 * ubiquitin_gpu.t.cc alone weren't enough -- see the "# totals" fixed-
 * width block below, whose field order must be kept in sync with
 * _print_energyred_helper by hand, same caveat that applied to the
 * abandoned Python version of this idea).
 */

#pragma once

#include <string>
#include <vector>

namespace check_ubiquitin {

  /**
   * Fixed-size "# totals" block layout inside every ENERGY03 block,
   * in on-disk order -- mirrors
   * io::Out_Configuration::_print_energyred_helper() exactly
   * (out_configuration.cc). Only the fields ubiquitin_runner.h's
   * StepResult already tracks are named; the rest are skipped over
   * by index, not reparsed into named fields, to avoid a second
   * 47-entry struct nobody would maintain in parallel.
   */
  constexpr int TOTALS_FIELD_COUNT = 47;
  constexpr int TOTALS_TOTAL = 0;
  constexpr int TOTALS_KINETIC_TOTAL = 1;
  constexpr int TOTALS_POTENTIAL_TOTAL = 2;
  constexpr int TOTALS_BONDED_TOTAL = 3;
  constexpr int TOTALS_NONBONDED_TOTAL = 9;
  constexpr int TOTALS_LJ_TOTAL = 10;
  constexpr int TOTALS_CRF_TOTAL = 11;
  constexpr int TOTALS_CONSTRAINTS_TOTAL = 23;

  struct TreBath {
    double total, com, ir, scale;
  };

  struct TreStep {
    unsigned step;
    double time;
    double totals[TOTALS_FIELD_COUNT];
    std::vector<TreBath> baths;
  };

  /**
   * @brief Parses every TIMESTEP/ENERGY03[/VOLUMEPRESSURE03] block
   * triple in `path`, in on-disk order. Returns false (and writes a
   * diagnostic to `err`) if the file can't be opened, contains fewer
   * than TOTALS_FIELD_COUNT numbers in some "# totals" block, or a
   * VOLUMEPRESSURE03 block's own structure doesn't parse -- a
   * malformed/truncated/garbled file is exactly the class of bug this
   * exists to catch, so a parse failure must be reported as loudly as
   * a value mismatch, never silently skipped.
   */
  bool parse_tre(const std::string & path, std::vector<TreStep> & out,
                  std::ostream & err);

  /**
   * @brief Scans every POSITIONRED block in a .trc trajectory file
   * (three whitespace-separated doubles per atom line) and checks each
   * one is finite and within a generously wide, physically-sane range
   * -- catches NaN/Inf and gross corruption (e.g. a formatting bug
   * that shifts columns) without pretending to be a real physics
   * check. Returns false (with a diagnostic to `err`) on the first
   * problem found, or if the file has no POSITIONRED blocks at all.
   */
  bool scan_trc_positions(const std::string & path, std::ostream & err);

}
