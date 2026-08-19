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
 * @file angle_gpu.t.cc
 * End-to-end correctness test for CUDA_Angle_Interaction (PLAN.md §10
 * step 14, bonded forces) -- runs the real CPU Angle_Interaction and the
 * real CUDA_Angle_Interaction on aladip's unperturbed topology/
 * configuration and compares the resulting per-atom forces,
 * per-energy-group angle energies, and atomic virial tensor. Only
 * USE_CUDA builds run this.
 */

#include "../stdheader.h"

#include <cmath>

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"

#include "../algorithm/algorithm/algorithm_sequence.h"
#include "../interaction/interaction.h"
#include "../interaction/bonded/angle_interaction.h"
#include "../interaction/bonded/cuda_angle_interaction.h"

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

    cpu_s.conf.current().force = 0.0;
    gpu_s.conf.current().force = 0.0;
    cpu_s.conf.current().virial_tensor = 0.0;
    gpu_s.conf.current().virial_tensor = 0.0;
    for (auto & e : cpu_s.conf.current().energies.angle_energy) e = 0.0;
    for (auto & e : gpu_s.conf.current().energies.angle_energy) e = 0.0;

    interaction::Angle_Interaction cpu_angle;
    interaction::CUDA_Angle_Interaction gpu_angle;

    if (cpu_angle.init(cpu_s.topo, cpu_s.conf, cpu_s.sim, std::cout, quiet) != 0 ||
        gpu_angle.init(gpu_s.topo, gpu_s.conf, gpu_s.sim, std::cout, quiet) != 0) {
      flush_messages("angle init");
      std::cerr << label << ": init() failed" << std::endl;
      return 1;
    }
    flush_messages("angle init");

    if (cpu_angle.calculate_interactions(cpu_s.topo, cpu_s.conf, cpu_s.sim) != 0 ||
        gpu_angle.calculate_interactions(gpu_s.topo, gpu_s.conf, gpu_s.sim) != 0) {
      flush_messages("angle calculate_interactions");
      std::cerr << label << ": calculate_interactions() failed" << std::endl;
      return 1;
    }
    flush_messages("angle calculate_interactions");

    // CUDA_Angle_Interaction writes force directly into the GPU mirror
    // and mark_gpu_dirty()s it instead of syncing back every call
    // (Algorithm_Sequence::run() does this publish automatically via
    // flush_gpu_dirty() before the next algorithm that needs it; this
    // standalone test calls calculate_interactions() directly, so it
    // must do the same publish itself before reading conf.current().
    // force() below).
    gpu_s.sim.cuda().flush_gpu_dirty(gpu_s.conf, gpu::MIRROR_FORCE);

    // See quartic_bond_gpu.t.cc's tolerance comment: bonded terms'
    // cost - cos0 (or dist2 - r0^2) formulas are near-cancellations near
    // equilibrium, amplifying the default FP_PRECISION=1 float position
    // mirror's truncation error well past nonbonded's 1e-4 tolerance.
    const double tol = 5e-3;
    int errors = 0;

    const unsigned num_atoms = static_cast<unsigned>(cpu_s.topo.num_atoms());
    for (unsigned i = 0; i < num_atoms; ++i) {
      const math::Vec diff = cpu_s.conf.current().force(i) - gpu_s.conf.current().force(i);
      const double scale = std::max(1.0, math::abs(cpu_s.conf.current().force(i)));
      if (math::abs(diff) > tol * scale) {
        std::cerr << label << ": force mismatch at atom " << i
                  << ": cpu=" << math::v2s(cpu_s.conf.current().force(i))
                  << " gpu=" << math::v2s(gpu_s.conf.current().force(i)) << std::endl;
        ++errors;
      }
    }

    const unsigned num_energy_groups =
        static_cast<unsigned>(cpu_s.conf.current().energies.angle_energy.size());
    for (unsigned g = 0; g < num_energy_groups; ++g) {
      const double cpu_e = cpu_s.conf.current().energies.angle_energy[g];
      const double gpu_e = gpu_s.conf.current().energies.angle_energy[g];
      if (std::abs(cpu_e - gpu_e) > tol * std::max(1.0, std::abs(cpu_e))) {
        std::cerr << label << ": angle_energy[" << g << "] mismatch: cpu=" << cpu_e
                  << " gpu=" << gpu_e << std::endl;
        ++errors;
      }
    }

    for (unsigned a = 0; a < 3; ++a) {
      for (unsigned b = 0; b < 3; ++b) {
        const double cpu_v = cpu_s.conf.current().virial_tensor(a, b);
        const double gpu_v = gpu_s.conf.current().virial_tensor(a, b);
        if (std::abs(cpu_v - gpu_v) > tol * std::max(1.0, std::abs(cpu_v))) {
          std::cerr << label << ": virial(" << a << "," << b << ") mismatch: cpu="
                    << cpu_v << " gpu=" << gpu_v << std::endl;
          ++errors;
        }
      }
    }

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

  return run_case(stopo, sconf, sinput, "angle_gpu", quiet);
}
