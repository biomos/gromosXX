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
 * @file ubiquitin_runner.h
 * Drives the real ubiquitin+water simulation (extended_test/, a full
 * ~22.7k-atom protein -- the largest, most structurally realistic
 * system exercised anywhere in this test suite) exactly the way
 * program/md does: io::read_input() + Algorithm_Sequence::init()/run(),
 * no shortcuts through a narrower test harness. Collects a handful of
 * key energies/bath kinetic energies directly from
 * configuration::Configuration and simulation::Multibath after every
 * step, in-memory -- no trajectory files are written or parsed (the
 * @fin/@trc/@tre arguments io::Out_Configuration would need are simply
 * never provided, since nothing here ever constructs one).
 */

#pragma once

#include <iostream>
#include <string>
#include <vector>

namespace check_ubiquitin {

  struct BathEnergy {
    double ekin;
  };

  struct StepResult {
    unsigned step;
    double total;
    double kinetic_total;
    double potential_total;
    double bonded_total;
    double nonbonded_total;
    double lj_total;
    double crf_total;
    double constraints_total;
    std::vector<BathEnergy> baths;
  };

  /**
   * @brief Runs `num_steps` real MD steps starting from the ubiquitin
   * topology/configuration with parameters from `imd_path`, collecting
   * one StepResult per step (including step 0, the initial force/
   * energy evaluation before any integration). Returns 0 on success;
   * on failure, writes a diagnostic to `err` and returns nonzero
   * (a hard error from io::read_input, Algorithm_Sequence::init(), or
   * any step's Algorithm_Sequence::run(), or a CRITICAL/ERROR-severity
   * io::messages entry logged along the way -- matching program/md's
   * own "Errors during initialisation!"/"Error during MD run!" checks,
   * so a silently-wrong run can't masquerade as success here either).
   */
  int run_ubiquitin(const std::string & imd_path, unsigned num_steps,
                     std::vector<StepResult> & out, std::ostream & err);

  /**
   * @brief Plain-text reference format, one line per step:
   * `step total kinetic_total potential_total bonded_total
   * nonbonded_total lj_total crf_total constraints_total num_baths
   * ekin[0] ekin[1] ...`. No JSON/serialization library needed --
   * this project doesn't otherwise depend on one, and the format is
   * only ever read back by write_reference()'s own counterpart below.
   */
  void write_reference(const std::string & path, const std::vector<StepResult> & steps);

  /**
   * @brief Returns false (and writes a diagnostic to `err`) if `path`
   * can't be opened or is malformed.
   */
  bool read_reference(const std::string & path, std::vector<StepResult> & out,
                       std::ostream & err);

}
