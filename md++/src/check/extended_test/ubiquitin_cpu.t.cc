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
 * @file ubiquitin_cpu.t.cc
 * CPU regression guard for the real ubiquitin+water simulation
 * (extended_test/, ~22.7k atoms -- the largest, most structurally
 * realistic system in this suite): re-runs the exact same 100 steps
 * ubiquitin_generate_reference.t.cc produced the committed reference
 * from, and requires every step's energies/bath kinetic energies to
 * match tightly. CPU is single-threaded and deterministic in this
 * build (no OpenMP linked), so this is effectively a bit-for-bit guard
 * against silent CPU-side regressions, not just a physics sanity
 * check. Registered unconditionally (not gated on USE_CUDA) -- nothing
 * here touches CUDA.
 */

#include "../../stdheader.h"
#include "ubiquitin_runner.h"

namespace {
  bool nearly_equal(double got, double ref, double rtol, double atol) {
    return std::abs(got - ref) <= atol + rtol * std::abs(ref);
  }
}

int main(int, char**) {
  std::vector<check_ubiquitin::StepResult> steps;
  const int rc = check_ubiquitin::run_ubiquitin(
      TOP_SOURCE_DIR "/src/check/extended_test/md_ubiquitin.imd", 100, steps, std::cerr);
  if (rc != 0) {
    std::cerr << "ubiquitin_cpu: run failed" << std::endl;
    return 1;
  }

  std::vector<check_ubiquitin::StepResult> reference;
  const std::string ref_path =
      TOP_SOURCE_DIR "/src/check/extended_test/reference/ubiquitin_cpu_reference.txt";
  if (!check_ubiquitin::read_reference(ref_path, reference, std::cerr)) {
    return 1;
  }

  if (steps.size() != reference.size()) {
    std::cerr << "ubiquitin_cpu: FAILED -- step count mismatch: got "
              << steps.size() << ", reference has " << reference.size() << std::endl;
    return 1;
  }

  const double rtol = 1e-6, atol = 1e-6;
  int errors = 0;
  for (size_t i = 0; i < steps.size(); ++i) {
    const check_ubiquitin::StepResult & g = steps[i];
    const check_ubiquitin::StepResult & r = reference[i];

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
      if (!nearly_equal(f.second.first, f.second.second, rtol, atol)) {
        std::cerr << "ubiquitin_cpu: step " << g.step << ": " << f.first
                  << " mismatch: got=" << f.second.first
                  << " ref=" << f.second.second << std::endl;
        ++errors;
      }
    }
    if (g.baths.size() != r.baths.size()) {
      std::cerr << "ubiquitin_cpu: step " << g.step << ": bath count mismatch: got="
                << g.baths.size() << " ref=" << r.baths.size() << std::endl;
      ++errors;
    } else {
      for (size_t b = 0; b < g.baths.size(); ++b) {
        if (!nearly_equal(g.baths[b].ekin, r.baths[b].ekin, rtol, atol)) {
          std::cerr << "ubiquitin_cpu: step " << g.step << ": bath " << b
                    << " ekin mismatch: got=" << g.baths[b].ekin
                    << " ref=" << r.baths[b].ekin << std::endl;
          ++errors;
        }
      }
    }
  }

  if (errors) {
    std::cerr << "ubiquitin_cpu: FAILED (" << errors << " mismatch(es))" << std::endl;
    return 1;
  }
  std::cout << "ubiquitin_cpu: OK (" << steps.size() << " step(s) match reference "
            << "within rtol=" << rtol << ")" << std::endl;
  return 0;
}
