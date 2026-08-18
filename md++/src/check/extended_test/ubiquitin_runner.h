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
 * step, in-memory.
 *
 * Also (optionally, via `out_paths`) drives a real
 * io::Out_Configuration exactly like program/md's own main loop
 * (traj.write() before each run(), traj.print() after, traj.write(...,
 * io::final)/traj.print_final() at the end) and writes real .fin/.trc/
 * .tre files -- this is deliberately a *second*, independent code path
 * from the in-memory collection above: a bug in Out_Configuration's
 * own formatting, precision, or which configuration object it reads
 * from would be invisible to a test that only inspects memory
 * directly. ubiquitin_cpu.t.cc/ubiquitin_gpu.t.cc cross-check the
 * parsed .tre contents (tre_parser.h) against these in-memory values
 * for exactly this reason.
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
   * @brief Real trajectory output file paths for run_ubiquitin() below
   * -- mirrors program/md's @fin/@trc/@tre arguments exactly (all
   * three are required together: io::Out_Configuration::init() reports
   * an error for any WRITETRAJ-requested trajectory with no matching
   * argument, same as a real run would).
   */
  struct OutputPaths {
    std::string fin;
    std::string trc;
    std::string tre;
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
   *
   * `out_paths`, if non-null, also drives a real io::Out_Configuration
   * writing actual .fin/.trc/.tre files at those paths (see this
   * file's own doc comment for why) -- the run is otherwise identical
   * either way (same steps, same in-memory collection into `out`).
   */
  int run_ubiquitin(const std::string & imd_path, unsigned num_steps,
                     std::vector<StepResult> & out, std::ostream & err,
                     const OutputPaths * out_paths = nullptr);

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

  /**
   * @brief Creates a fresh, empty, process-unique scratch directory
   * (under the system temp dir) and returns its path, or "" (with a
   * diagnostic to `err`) on failure. Callers are responsible for
   * cleaning it up (or not -- it's under the OS temp dir either way).
   */
  std::string make_temp_dir(std::ostream & err);

}
