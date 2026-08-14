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
 * @file leap_frog_gpu.t.cc
 * Validates Leap_Frog_Velocity<gpuBackend>/Leap_Frog_Position<gpuBackend>
 * (PLAN.md §10 roadmap step 11) against the CPU reference implementation,
 * driven over several synthetic steps with deterministic, externally
 * supplied forces (no real nonbonded evaluation needed -- this test is
 * about the integrator, not the force field). Two independently-built
 * simulations, identical topology/config, fed the identical per-step
 * force pattern; positions/velocities must match within float tolerance
 * (the GPU path runs in FPL_TYPE, typically single precision) at every
 * step.
 *
 * Only USE_CUDA builds run this.
 */

#include "../stdheader.h"

#include <cmath>

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"

#include "../algorithm/algorithm/algorithm_sequence.h"
#include "../algorithm/integration/leap_frog.h"

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

  void set_deterministic_force(configuration::Configuration & conf, unsigned num_atoms,
                                unsigned step) {
    for (unsigned i = 0; i < num_atoms; ++i) {
      const double a = static_cast<double>(i) + static_cast<double>(step) * 0.37;
      conf.current().force(i) = math::Vec(std::sin(a * 0.7 + 1.0),
                                           std::cos(a * 1.3 + 2.0),
                                           std::sin(a * 2.1 + 3.0));
    }
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

    const unsigned num_atoms = static_cast<unsigned>(cpu_s.topo.num_atoms());
    const unsigned num_steps = 5;

    algorithm::Leap_Frog_Velocity<util::cpuBackend> cpu_velocity;
    algorithm::Leap_Frog_Position<util::cpuBackend> cpu_position;
    algorithm::Leap_Frog_Velocity<util::gpuBackend> gpu_velocity;
    algorithm::Leap_Frog_Position<util::gpuBackend> gpu_position;

    int errors = 0;
    const double tol = 1e-4;

    for (unsigned step = 0; step < num_steps; ++step) {
      cpu_s.sim.steps() = step;
      gpu_s.sim.steps() = step;

      set_deterministic_force(cpu_s.conf, num_atoms, step);
      set_deterministic_force(gpu_s.conf, num_atoms, step);

      if (cpu_velocity.apply(cpu_s.topo, cpu_s.conf, cpu_s.sim) != 0 ||
          cpu_position.apply(cpu_s.topo, cpu_s.conf, cpu_s.sim) != 0) {
        std::cerr << label << ": CPU leap-frog apply() failed at step " << step << std::endl;
        ++errors;
        break;
      }
      if (gpu_velocity.apply(gpu_s.topo, gpu_s.conf, gpu_s.sim) != 0 ||
          gpu_position.apply(gpu_s.topo, gpu_s.conf, gpu_s.sim) != 0) {
        flush_messages("gpu leap-frog apply");
        std::cerr << label << ": GPU leap-frog apply() failed at step " << step << std::endl;
        ++errors;
        break;
      }
      flush_messages("gpu leap-frog apply");

      for (unsigned i = 0; i < num_atoms; ++i) {
        const math::Vec dv = cpu_s.conf.current().vel(i) - gpu_s.conf.current().vel(i);
        const math::Vec dx = cpu_s.conf.current().pos(i) - gpu_s.conf.current().pos(i);
        const double vscale = std::max(1.0, math::abs(cpu_s.conf.current().vel(i)));
        const double xscale = std::max(1.0, math::abs(cpu_s.conf.current().pos(i)));
        if (math::abs(dv) > tol * vscale) {
          std::cerr << label << ": vel mismatch at step " << step << ", atom " << i
                    << ": cpu=" << math::v2s(cpu_s.conf.current().vel(i))
                    << " gpu=" << math::v2s(gpu_s.conf.current().vel(i)) << std::endl;
          ++errors;
        }
        if (math::abs(dx) > tol * xscale) {
          std::cerr << label << ": pos mismatch at step " << step << ", atom " << i
                    << ": cpu=" << math::v2s(cpu_s.conf.current().pos(i))
                    << " gpu=" << math::v2s(gpu_s.conf.current().pos(i)) << std::endl;
          ++errors;
        }
      }
    }

    if (errors == 0) {
      std::cout << label << ": OK (" << num_steps << " steps, " << num_atoms
                << " atoms, CPU/GPU leap-frog match)" << std::endl;
    } else {
      std::cerr << label << ": FAILED (" << errors << " mismatch(es))" << std::endl;
    }
    return errors == 0 ? 0 : 1;
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

  return run_case(stopo, sconf, sinput, "leap_frog_gpu", quiet);
}
