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
 * @file remove_com_motion_gpu.t.cc
 * Validates Remove_COM_Motion<gpuBackend> (PLAN.md §10 step 13 --
 * broadening algorithm coverage) against the CPU reference
 * (remove_com_motion_cpu.cc): translation removal, rotation removal, and
 * the combined apply() path, on aladip's real topology/config. Only
 * USE_CUDA builds run this.
 */

#include "../stdheader.h"

#include <cmath>

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"

#include "../algorithm/algorithm/algorithm_sequence.h"
#include "../algorithm/constraints/remove_com_motion.h"

#include "../io/argument.h"
#include "../util/parse_verbosity.h"
#include "../util/usage.h"
#include "../io/topology/in_topology.h"
#include "../io/message.h"

#include "../util/create_simulation.h"

#include "check.h"

#ifdef XXMPI
  #include <mpi.h>
#endif

namespace {

  bool flush_messages(const char * label) {
    const io::message::severity_enum sev = io::messages.display(std::cerr);
    io::messages.clear();
    return sev >= io::message::error;
  }

  int compare_velocities(const util::simulation_struct & cpu_s,
                          const util::simulation_struct & gpu_s,
                          const char * label, double tol) {
    int errors = 0;
    const unsigned num_atoms = static_cast<unsigned>(cpu_s.topo.num_atoms());
    for (unsigned i = 0; i < num_atoms; ++i) {
      const math::Vec diff = cpu_s.conf.current().vel(i) - gpu_s.conf.current().vel(i);
      const double scale = std::max(1.0, math::abs(cpu_s.conf.current().vel(i)));
      if (math::abs(diff) > tol * scale) {
        std::cerr << label << ": vel mismatch at atom " << i
                  << ": cpu=" << math::v2s(cpu_s.conf.current().vel(i))
                  << " gpu=" << math::v2s(gpu_s.conf.current().vel(i)) << std::endl;
        ++errors;
      }
    }
    return errors;
  }

  int run_case(const std::string & stopo, const std::string & sconf,
               const std::string & sinput, const char * label, bool quiet) {
    util::simulation_struct cpu_s, gpu_s;
    io::In_Topology cpu_in_topo, gpu_in_topo;
    cpu_in_topo.quiet = quiet;
    gpu_in_topo.quiet = quiet;

    if (util::create_simulation(stopo, "", sconf, sinput, cpu_s, cpu_in_topo,
                                 "", "", "", "", "", "", "", "", quiet) != 0 ||
        util::create_simulation(stopo, "", sconf, sinput, gpu_s, gpu_in_topo,
                                 "", "", "", "", "", "", "", "", quiet) != 0) {
      std::cerr << label << ": creating simulation failed" << std::endl;
      return 1;
    }
    io::messages.display(std::cout);
    io::messages.clear();

    cpu_s.sim.param().centreofmass.remove_trans = true;
    cpu_s.sim.param().centreofmass.remove_rot   = true;
    cpu_s.sim.steps() = 1;
    gpu_s.sim.param().centreofmass.remove_trans = true;
    gpu_s.sim.param().centreofmass.remove_rot   = true;
    gpu_s.sim.steps() = 1;

    const double tol = 1e-4;
    int errors = 0;

    // --- direct comparison of remove_com_translation()'s return value
    // and velocity effect ---
    algorithm::Remove_COM_Motion<util::cpuBackend> cpu_rcom;
    algorithm::Remove_COM_Motion<util::gpuBackend> gpu_rcom;

    const double cpu_ekin_trans = cpu_rcom.remove_com_translation(cpu_s.topo, cpu_s.conf, cpu_s.sim, true);
    const double gpu_ekin_trans = gpu_rcom.remove_com_translation(gpu_s.topo, gpu_s.conf, gpu_s.sim, true);
    flush_messages("gpu remove_com_translation");
    // remove_com_translation() leaves VEL dirty-but-unpublished on the
    // GPU mirror (mark_gpu_dirty(), no CPU round trip) -- this test
    // reads conf.current().vel directly, so it must ask for the
    // publish itself, same as leap_frog_gpu.t.cc's identical fix.
    gpu_s.sim.cuda().flush_gpu_dirty(gpu_s.conf, gpu::MIRROR_VEL);
    if (std::abs(cpu_ekin_trans - gpu_ekin_trans) > tol * std::max(1.0, std::abs(cpu_ekin_trans))) {
      std::cerr << label << ": ekin_trans mismatch: cpu=" << cpu_ekin_trans
                << " gpu=" << gpu_ekin_trans << std::endl;
      ++errors;
    }
    errors += compare_velocities(cpu_s, gpu_s, (std::string(label) + " (translation)").c_str(), tol);

    // --- rotation removal, continuing from the post-translation state ---
    const double cpu_ekin_rot = cpu_rcom.remove_com_rotation(cpu_s.topo, cpu_s.conf, cpu_s.sim, true);
    const double gpu_ekin_rot = gpu_rcom.remove_com_rotation(gpu_s.topo, gpu_s.conf, gpu_s.sim, true);
    flush_messages("gpu remove_com_rotation");
    gpu_s.sim.cuda().flush_gpu_dirty(gpu_s.conf, gpu::MIRROR_POS | gpu::MIRROR_VEL);
    if (std::abs(cpu_ekin_rot - gpu_ekin_rot) > tol * std::max(1.0, std::abs(cpu_ekin_rot))) {
      std::cerr << label << ": ekin_rot mismatch: cpu=" << cpu_ekin_rot
                << " gpu=" << gpu_ekin_rot << std::endl;
      ++errors;
    }
    errors += compare_velocities(cpu_s, gpu_s, (std::string(label) + " (rotation)").c_str(), tol);

    if (errors) {
      std::cerr << label << ": FAILED (" << errors << " mismatch(es))" << std::endl;
    } else {
      std::cout << label << ": OK (ekin_trans=" << gpu_ekin_trans
                << " ekin_rot=" << gpu_ekin_rot << ")" << std::endl;
    }
    return errors ? 1 : 0;
  }

  /**
   * Drives the real Algorithm::apply() path (skip_step-gated, matching
   * how create_md_sequence.cc actually calls this every step) rather
   * than the two removal methods directly -- confirms the full
   * dispatch, not just the underlying math.
   */
  int run_apply_case(const std::string & stopo, const std::string & sconf,
                      const std::string & sinput, const char * label, bool quiet) {
    util::simulation_struct cpu_s, gpu_s;
    io::In_Topology cpu_in_topo, gpu_in_topo;
    cpu_in_topo.quiet = quiet;
    gpu_in_topo.quiet = quiet;

    if (util::create_simulation(stopo, "", sconf, sinput, cpu_s, cpu_in_topo,
                                 "", "", "", "", "", "", "", "", quiet) != 0 ||
        util::create_simulation(stopo, "", sconf, sinput, gpu_s, gpu_in_topo,
                                 "", "", "", "", "", "", "", "", quiet) != 0) {
      std::cerr << label << ": creating simulation failed" << std::endl;
      return 1;
    }
    io::messages.display(std::cout);
    io::messages.clear();

    for (util::simulation_struct * s : {&cpu_s, &gpu_s}) {
      s->sim.param().centreofmass.skip_step = 1;
      s->sim.param().centreofmass.remove_trans = true;
      s->sim.param().centreofmass.remove_rot   = true;
      s->sim.steps() = 1;
    }

    algorithm::Remove_COM_Motion<util::cpuBackend> cpu_rcom;
    algorithm::Remove_COM_Motion<util::gpuBackend> gpu_rcom;

    if (cpu_rcom.apply(cpu_s.topo, cpu_s.conf, cpu_s.sim) != 0 ||
        gpu_rcom.apply(gpu_s.topo, gpu_s.conf, gpu_s.sim) != 0) {
      flush_messages("apply");
      std::cerr << label << ": apply() failed" << std::endl;
      return 1;
    }
    flush_messages("apply");
    // Same as above -- apply() defers the publish; ask for it before
    // reading conf.current().vel directly.
    gpu_s.sim.cuda().flush_gpu_dirty(gpu_s.conf, gpu::MIRROR_POS | gpu::MIRROR_VEL);

    const int errors = compare_velocities(cpu_s, gpu_s, label, 1e-4);
    if (errors) {
      std::cerr << label << ": FAILED (" << errors << " mismatch(es))" << std::endl;
    } else {
      std::cout << label << ": OK" << std::endl;
    }
    return errors ? 1 : 0;
  }

} // namespace

int main(int argc, char* argv[]) {
#ifdef XXMPI
  MPI_Init(&argc, &argv);
#endif

  util::Known knowns;
  knowns << "verb";

  std::string usage = argv[0];
  usage += "\n\t[@verb   <[module:][submodule:]level>]\n";

  io::Argument args;
  if (args.parse(argc, argv, knowns, true)) {
    std::cerr << usage << std::endl;
    return 1;
  }
  util::parse_verbosity(args);
  const bool quiet = (args.count("verb") == -1);

  std::string stopo, sconf, sinput;
  GETFILEPATH(stopo, "aladip.topo", "src/check/data/");
  GETFILEPATH(sconf, "aladip.conf", "src/check/data/");
  GETFILEPATH(sinput, "aladip_unperturbed.in", "src/check/data/");

  int total = 0;
  total += run_case(stopo, sconf, sinput, "remove_com_motion_gpu", quiet);
  total += run_apply_case(stopo, sconf, sinput, "remove_com_motion_gpu_apply", quiet);
  return total;
}
