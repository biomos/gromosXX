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
 * @file ubiquitin_gpu.t.cc
 * GPU counterpart of ubiquitin_cpu.t.cc, for the same real ubiquitin+
 * water simulation (~22.7k atoms) -- but a GPU MD trajectory is a
 * chaotic nonlinear dynamical system's realization under different
 * floating-point summation order (atomics) and mixed precision than
 * CPU's; it is *expected* to diverge from the CPU trajectory well
 * before 100 steps, no matter how correct the GPU code is. This test
 * therefore checks several different things, each appropriately:
 *
 *  - Step 0 (the initial force/energy evaluation, before any
 *    integration): CPU and GPU see the *identical* starting
 *    configuration here, so this is a real, direct correctness check
 *    of the GPU force field against the CPU reference -- not yet
 *    chaos-diverged. Checked with a moderate (not bit-for-bit)
 *    tolerance, since even this one comparison involves mixed-
 *    precision (float) accumulation over thousands of pairwise terms
 *    that CPU sums in double. Checked against both the in-memory
 *    energies *and* the actual .tre file contents (see below).
 *  - The rest of the trajectory: only a physical-sanity check (no
 *    NaN/Inf, every bath's kinetic energy stays in a believable,
 *    bounded range) -- not compared against the CPU trajectory at all.
 *  - The actual .tre/.trc file contents (io::Out_Configuration, driven
 *    exactly like program/md's own main loop) are cross-checked
 *    against the in-memory values above, and .trc positions are
 *    scanned for NaN/Inf/gross corruption -- catches a bug in the
 *    file-writing path itself, which is invisible to a test that only
 *    inspects memory directly.
 *
 * Only registered under USE_CUDA (see CMakeLists.txt) -- there is no
 * GPU accelerator to select otherwise.
 */

#include "../../stdheader.h"
#include "ubiquitin_runner.h"
#include "tre_parser.h"

namespace {
  bool nearly_equal(double got, double ref, double rtol, double atol) {
    return std::abs(got - ref) <= atol + rtol * std::abs(ref);
  }
}

int main(int, char**) {
  const std::string workdir = check_ubiquitin::make_temp_dir(std::cerr);
  if (workdir.empty()) return 1;
  check_ubiquitin::OutputPaths out_paths;
  out_paths.fin = workdir + "/ubiquitin.fin.cnf";
  out_paths.trc = workdir + "/ubiquitin.trc";
  out_paths.tre = workdir + "/ubiquitin.tre";

  std::vector<check_ubiquitin::StepResult> steps;
  const int rc = check_ubiquitin::run_ubiquitin(
      TOP_SOURCE_DIR "/src/check/extended_test/md_ubiquitin_gpu.imd", 100, steps,
      std::cerr, &out_paths);
  if (rc != 0) {
    std::cerr << "ubiquitin_gpu: run failed" << std::endl;
    return 1;
  }
  if (steps.empty()) {
    std::cerr << "ubiquitin_gpu: FAILED -- no steps produced" << std::endl;
    return 1;
  }

  std::vector<check_ubiquitin::StepResult> reference;
  const std::string ref_path =
      TOP_SOURCE_DIR "/src/check/extended_test/reference/ubiquitin_cpu_reference.txt";
  if (!check_ubiquitin::read_reference(ref_path, reference, std::cerr)) {
    return 1;
  }
  if (reference.empty()) {
    std::cerr << "ubiquitin_gpu: FAILED -- empty reference" << std::endl;
    return 1;
  }

  int errors = 0;

  // Step 0: real correctness check against the CPU reference.
  {
    const double rtol = 5e-3, atol = 1.0;
    const check_ubiquitin::StepResult & g = steps[0];
    const check_ubiquitin::StepResult & r = reference[0];
    const std::pair<const char *, std::pair<double, double> > fields[] = {
      {"total",            {g.total,            r.total}},
      {"kinetic_total",    {g.kinetic_total,    r.kinetic_total}},
      {"potential_total",  {g.potential_total,  r.potential_total}},
      {"bonded_total",     {g.bonded_total,     r.bonded_total}},
      {"nonbonded_total",  {g.nonbonded_total,  r.nonbonded_total}},
      {"lj_total",         {g.lj_total,         r.lj_total}},
      {"crf_total",        {g.crf_total,        r.crf_total}},
      {"constraints_total",{g.constraints_total,r.constraints_total}},
    };
    for (const auto & f : fields) {
      const double diff = std::abs(f.second.first - f.second.second);
      if (diff > atol + rtol * std::abs(f.second.second)) {
        std::cerr << "ubiquitin_gpu: step 0: " << f.first
                  << " mismatch: got=" << f.second.first
                  << " ref=" << f.second.second << std::endl;
        ++errors;
      }
    }
  }

  // Whole trajectory: physical sanity, not numerical agreement with CPU.
  for (const check_ubiquitin::StepResult & s : steps) {
    const double vals[] = {s.total, s.kinetic_total, s.potential_total,
                            s.bonded_total, s.nonbonded_total, s.lj_total,
                            s.crf_total, s.constraints_total};
    for (double v : vals) {
      if (std::isnan(v) || std::isinf(v)) {
        std::cerr << "ubiquitin_gpu: step " << s.step << ": NaN/Inf energy ("
                  << v << ")" << std::endl;
        ++errors;
      }
    }
    for (const check_ubiquitin::BathEnergy & b : s.baths) {
      // A believable, generous band -- not a tight physical bound: the
      // point is to catch a blown-up/diverged/NaN-adjacent simulation,
      // not to validate thermostat accuracy (a separate, existing test,
      // nosehoover_gpu.t.cc/temperature_gpu.t.cc, already does that).
      if (std::isnan(b.ekin) || std::isinf(b.ekin) || b.ekin < 0.0 || b.ekin > 1.0e7) {
        std::cerr << "ubiquitin_gpu: step " << s.step
                  << ": bath kinetic energy out of sane range (" << b.ekin << ")"
                  << std::endl;
        ++errors;
      }
    }
  }

  if (errors) {
    std::cerr << "ubiquitin_gpu: FAILED (" << errors << " issue(s))" << std::endl;
    return 1;
  }

  // The actual .tre file, produced by the same real
  // io::Out_Configuration path program/md uses -- must agree with the
  // in-memory values above at every step that got written, and step 0
  // (the one point directly comparable to CPU) must still be within
  // the same moderate tolerance against the CPU reference when read
  // back from disk, not just in memory.
  std::vector<check_ubiquitin::TreStep> tre_steps;
  if (!check_ubiquitin::parse_tre(out_paths.tre, tre_steps, std::cerr)) {
    std::cerr << "ubiquitin_gpu: FAILED -- could not parse " << out_paths.tre
              << std::endl;
    return 1;
  }

  const double file_rtol = 1e-6, file_atol = 1e-3;
  int file_errors = 0;
  for (const check_ubiquitin::TreStep & ts : tre_steps) {
    if (ts.step >= steps.size()) {
      std::cerr << "ubiquitin_gpu: .tre has step " << ts.step
                << " beyond the " << steps.size() << " steps run" << std::endl;
      ++file_errors;
      continue;
    }
    const check_ubiquitin::StepResult & mem = steps[ts.step];
    const std::pair<const char *, std::pair<double, double> > fields[] = {
      {"total",            {ts.totals[check_ubiquitin::TOTALS_TOTAL],            mem.total}},
      {"kinetic_total",    {ts.totals[check_ubiquitin::TOTALS_KINETIC_TOTAL],    mem.kinetic_total}},
      {"potential_total",  {ts.totals[check_ubiquitin::TOTALS_POTENTIAL_TOTAL],  mem.potential_total}},
      {"bonded_total",     {ts.totals[check_ubiquitin::TOTALS_BONDED_TOTAL],     mem.bonded_total}},
      {"nonbonded_total",  {ts.totals[check_ubiquitin::TOTALS_NONBONDED_TOTAL],  mem.nonbonded_total}},
      {"lj_total",         {ts.totals[check_ubiquitin::TOTALS_LJ_TOTAL],         mem.lj_total}},
      {"crf_total",        {ts.totals[check_ubiquitin::TOTALS_CRF_TOTAL],        mem.crf_total}},
      {"constraints_total",{ts.totals[check_ubiquitin::TOTALS_CONSTRAINTS_TOTAL],mem.constraints_total}},
    };
    for (const auto & f : fields) {
      if (!nearly_equal(f.second.first, f.second.second, file_rtol, file_atol)) {
        std::cerr << "ubiquitin_gpu: .tre step " << ts.step << ": " << f.first
                  << " disagrees with in-memory value: file=" << f.second.first
                  << " memory=" << f.second.second << std::endl;
        ++file_errors;
      }
    }
  }

  if (!check_ubiquitin::scan_trc_positions(out_paths.trc, std::cerr)) {
    std::cerr << "ubiquitin_gpu: FAILED -- .trc position sanity check failed"
              << std::endl;
    return 1;
  }

  if (file_errors) {
    std::cerr << "ubiquitin_gpu: FAILED (" << file_errors
              << " .tre/in-memory mismatch(es))" << std::endl;
    return 1;
  }

  std::cout << "ubiquitin_gpu: OK (step 0 matches CPU reference, " << steps.size()
            << " step(s) stable, no NaN/Inf, " << tre_steps.size()
            << " .tre step(s) match in-memory values, .trc positions sane)"
            << std::endl;
  return 0;
}
