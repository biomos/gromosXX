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
 * @file ubiquitin_generate_reference.t.cc
 * NOT a ctest test (not registered via add_test) -- a manual, one-off
 * tool that runs the real ubiquitin+water simulation on CPU for 100
 * steps and writes extended_test/reference/ubiquitin_cpu_reference.txt,
 * the ground truth ubiquitin_cpu.t.cc/ubiquitin_gpu.t.cc compare
 * against. Re-run deliberately (and review the resulting file's diff
 * before committing it) whenever the underlying physics genuinely
 * changes -- a real bug fix, or an edit to md_ubiquitin.imd/
 * ubiquitin_54a8.top/ubiquitin.cnf.
 *
 * CPU is single-threaded and deterministic here (no OpenMP linked in
 * this build), so this reference is expected to reproduce bit-for-bit
 * across runs on the same machine/toolchain -- ubiquitin_cpu.t.cc
 * checks exactly that, with a tight tolerance.
 */

#include "../../stdheader.h"
#include "ubiquitin_runner.h"

int main(int, char**) {
  std::vector<check_ubiquitin::StepResult> steps;
  const int rc = check_ubiquitin::run_ubiquitin(
      TOP_SOURCE_DIR "/src/check/extended_test/md_ubiquitin.imd", 100, steps, std::cerr);
  if (rc != 0) {
    std::cerr << "ubiquitin_generate_reference: run failed" << std::endl;
    return 1;
  }
  const std::string out_path =
      TOP_SOURCE_DIR "/src/check/extended_test/reference/ubiquitin_cpu_reference.txt";
  check_ubiquitin::write_reference(out_path, steps);
  std::cout << "wrote " << out_path << " (" << steps.size() << " step(s))" << std::endl;
  return 0;
}
