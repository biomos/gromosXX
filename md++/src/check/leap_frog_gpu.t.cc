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
      // Stand-in for Forcefield::apply(): in a real sequence, the
      // algorithm that just (re-)computed force is followed by
      // Algorithm_Sequence::run()'s default post-apply() invalidation
      // (gpu_mirror_touches() == MIRROR_ALL), which is what tells the
      // GPU-mirror freshness tracker this step's force is new and must
      // be resynced before Leap_Frog_Velocity<gpuBackend> reads it.
      // This test drives Velocity/Position directly (no real
      // Algorithm_Sequence), so it must supply that invalidation itself.
      gpu_s.sim.cuda().invalidate_gpu_mirror(gpu_s.conf);

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

  /**
   * Stand-in for Berendsen_Thermostat/NoseHoover_Thermostat: a CPU-only
   * algorithm that rescales conf.current().vel directly, with no GPU
   * backend and no gpu_mirror_touches() override (uses the base class
   * default, MIRROR_ALL). This is exactly the shape of algorithm
   * create_md_sequence.cc inserts between Leap_Frog_Velocity and
   * Leap_Frog_Position whenever sim.param().multibath.couple is set --
   * the scenario that silently broke Leap_Frog_Position<gpuBackend>
   * before the GPU-mirror freshness tracking fix (it would have kept
   * reading the pre-scaling velocity straight off the GPU mirror).
   */
  class Vel_Scale_Algorithm : public algorithm::Algorithm {
  public:
    explicit Vel_Scale_Algorithm(double factor)
      : algorithm::Algorithm("Vel_Scale"), m_factor(factor) {}

    int init(topology::Topology &, configuration::Configuration &,
             simulation::Simulation &, std::ostream &, bool) override {
      return 0;
    }

    int apply(topology::Topology & topo, configuration::Configuration & conf,
               simulation::Simulation &) override {
      const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
      for (unsigned i = 0; i < num_atoms; ++i)
        conf.current().vel(i) *= m_factor;
      return 0;
    }

  private:
    double m_factor;
  };

  /**
   * Regression test for the bug this design fixes: drives
   * Leap_Frog_Velocity<gpuBackend> -> Vel_Scale_Algorithm (simulated
   * thermostat) -> Leap_Frog_Position<gpuBackend> through a real
   * algorithm::Algorithm_Sequence::run(), so the actual
   * invalidate_gpu_mirror() hook fires after each step, exactly as it
   * would in create_md_sequence.cc's real sequence. Before the fix,
   * Leap_Frog_Position<gpuBackend> would have silently used the
   * pre-scaling velocity (never invalidated, never resynced); after,
   * it must match a CPU-only run of the identical three-algorithm
   * sequence exactly.
   */
  int run_thermostat_interleaved_case(const std::string & stopo, const std::string & sconf,
                                       const std::string & sinput, const char * label,
                                       bool quiet) {
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
    const unsigned num_steps = 3;
    const double scale_factor = 0.97; // arbitrary, thermostat-like rescale

    algorithm::Algorithm_Sequence cpu_seq;
    cpu_seq.push_back(new algorithm::Leap_Frog_Velocity<util::cpuBackend>());
    cpu_seq.push_back(new Vel_Scale_Algorithm(scale_factor));
    cpu_seq.push_back(new algorithm::Leap_Frog_Position<util::cpuBackend>());

    algorithm::Algorithm_Sequence gpu_seq;
    gpu_seq.push_back(new algorithm::Leap_Frog_Velocity<util::gpuBackend>());
    gpu_seq.push_back(new Vel_Scale_Algorithm(scale_factor));
    gpu_seq.push_back(new algorithm::Leap_Frog_Position<util::gpuBackend>());

    int errors = 0;
    const double tol = 1e-4;

    for (unsigned step = 0; step < num_steps; ++step) {
      cpu_s.sim.steps() = step;
      gpu_s.sim.steps() = step;

      set_deterministic_force(cpu_s.conf, num_atoms, step);
      set_deterministic_force(gpu_s.conf, num_atoms, step);
      // Stand-in for Forcefield::apply()'s own post-apply()
      // invalidation within a real Algorithm_Sequence -- see the
      // identical comment in run_case() above.
      gpu_s.sim.cuda().invalidate_gpu_mirror(gpu_s.conf);

      if (cpu_seq.run(cpu_s.topo, cpu_s.conf, cpu_s.sim) != 0) {
        std::cerr << label << ": CPU sequence failed at step " << step << std::endl;
        ++errors;
        break;
      }
      if (gpu_seq.run(gpu_s.topo, gpu_s.conf, gpu_s.sim) != 0) {
        flush_messages("gpu thermostat-interleaved sequence");
        std::cerr << label << ": GPU sequence failed at step " << step << std::endl;
        ++errors;
        break;
      }
      flush_messages("gpu thermostat-interleaved sequence");

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
                << " atoms, thermostat-interleaved CPU/GPU sequence match)" << std::endl;
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

  int total = 0;
  total += run_case(stopo, sconf, sinput, "leap_frog_gpu", quiet);
  total += run_thermostat_interleaved_case(stopo, sconf, sinput,
                                            "leap_frog_gpu_thermostat_interleaved", quiet);
  return total;
}
