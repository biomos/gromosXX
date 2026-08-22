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
 * @file settle_gpu.t.cc
 * End-to-end correctness test for CUDA_Settle (PLAN.md §10 step 18,
 * constraints) -- runs the real CPU algorithm::Settle and the real
 * algorithm::CUDA_Settle on aladip's unperturbed topology/
 * configuration (its solvent is 3-site SPC water, exactly CUDA_Settle's
 * v1 scope) and compares the resulting positions, velocities,
 * constraint forces, and virial tensor. Same displaced-solvent-
 * positions scenario as shake_gpu.t.cc: aladip's own starting
 * configuration already nearly satisfies every constraint. Only
 * USE_CUDA builds run this.
 */

#include "../stdheader.h"

#include <cmath>

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"

#include "../algorithm/algorithm/algorithm_sequence.h"
#include "../algorithm/constraints/settle.h"
#include "../algorithm/constraints/cuda_settle.h"
#include "gpu/cuda/manager/cuda_manager.h"

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

  void displace_solvent(util::simulation_struct & s) {
    const unsigned num_atoms = static_cast<unsigned>(s.topo.num_atoms());
    for (unsigned i = static_cast<unsigned>(s.topo.num_solute_atoms()); i < num_atoms; ++i) {
      const double d = 0.005 * std::sin(0.7 * i + 1.0);
      s.conf.current().pos(i) += math::Vec(d, -0.5 * d, 0.25 * d);
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

    // Drive SETTLE regardless of the input file's own NTCS choice
    // (aladip_unperturbed.in uses SHAKE for solvent -- see
    // shake_gpu.t.cc); system.nsm must be set for Settle::apply() to
    // actually run (it gates on `sim.param().system.nsm &&
    // constraint.solvent.algorithm == constr_settle`).
    cpu_s.sim.param().constraint.solvent.algorithm = simulation::constr_settle;
    gpu_s.sim.param().constraint.solvent.algorithm = simulation::constr_settle;
    cpu_s.sim.param().system.nsm = static_cast<int>(cpu_s.topo.num_solvent_molecules(0));
    gpu_s.sim.param().system.nsm = static_cast<int>(gpu_s.topo.num_solvent_molecules(0));

    cpu_s.conf.old() = cpu_s.conf.current();
    gpu_s.conf.old() = gpu_s.conf.current();

    displace_solvent(cpu_s);
    displace_solvent(gpu_s);

    algorithm::Settle cpu_settle;
    algorithm::CUDA_Settle gpu_settle;

    if (cpu_settle.init(cpu_s.topo, cpu_s.conf, cpu_s.sim, std::cout, quiet) != 0 ||
        gpu_settle.init(gpu_s.topo, gpu_s.conf, gpu_s.sim, std::cout, quiet) != 0) {
      flush_messages("settle init");
      std::cerr << label << ": init() failed" << std::endl;
      return 1;
    }
    flush_messages("settle init");

    const int cpu_rc = cpu_settle.apply(cpu_s.topo, cpu_s.conf, cpu_s.sim);
    const int gpu_rc = gpu_settle.apply(gpu_s.topo, gpu_s.conf, gpu_s.sim);
    flush_messages("settle apply");

    // CUDA_Settle leaves pos/vel GPU-resident (mark_gpu_dirty(), no
    // eager CPU publish) -- this test reads gpu_s.conf directly, so it
    // must request the publish explicitly, same as any other direct
    // consumer.
    gpu_s.sim.cuda().flush_gpu_dirty(gpu_s.conf, gpu::MIRROR_POS | gpu::MIRROR_VEL);

    if (cpu_rc != 0 || gpu_rc != 0) {
      std::cerr << label << ": apply() failed (cpu_rc=" << cpu_rc
                << " gpu_rc=" << gpu_rc << ")" << std::endl;
      return 1;
    }

    // SETTLE is a single closed-form evaluation per molecule (no
    // iteration) -- both sides run the exact same arithmetic sequence,
    // so pos itself stays tight (bit-comparable up to associativity).
    // vel/constraint_force need a looser, scale-relative tolerance
    // (same shape as shake_gpu.t.cc's, for the same reason): CUDA_
    // Settle now reads/writes pos through the shared GPU mirror
    // (FPL/float precision, matching every other GPU-resident
    // algorithm) instead of a private full-double upload/download --
    // settle_kernels.cu's own math still runs in double, but its input
    // position already carries the mirror's float rounding. SETTLE's
    // correction delta (`d_a` in settle_kernels.cu) is a *difference*
    // of two close double3 values derived from that float-rounded
    // input -- ordinary catastrophic cancellation amplifies that
    // rounding into a much larger relative error in the delta itself,
    // and vel/constraint_force are both directly proportional to that
    // delta (scaled by 1/dt and 1/dt^2 respectively).
    const double tol = 1e-6;
    const double cf_atol = 1.0, cf_rtol = 2e-3;
    const double vel_atol = 2e-3, vel_rtol = 5e-4;
    int errors = 0;

    const unsigned num_atoms = static_cast<unsigned>(cpu_s.topo.num_atoms());
    const unsigned first_solvent = static_cast<unsigned>(cpu_s.topo.num_solute_atoms());
    for (unsigned i = first_solvent; i < num_atoms; ++i) {
      const math::Vec pos_diff = cpu_s.conf.current().pos(i) - gpu_s.conf.current().pos(i);
      const double pos_scale = std::max(1.0, math::abs(cpu_s.conf.current().pos(i)));
      if (math::abs(pos_diff) > tol * pos_scale) {
        std::cerr << label << ": pos mismatch at atom " << i
                  << ": cpu=" << math::v2s(cpu_s.conf.current().pos(i))
                  << " gpu=" << math::v2s(gpu_s.conf.current().pos(i)) << std::endl;
        ++errors;
      }
      const math::Vec vel_diff = cpu_s.conf.current().vel(i) - gpu_s.conf.current().vel(i);
      if (math::abs(vel_diff) > vel_atol + vel_rtol * math::abs(cpu_s.conf.current().vel(i))) {
        std::cerr << label << ": vel mismatch at atom " << i
                  << ": cpu=" << math::v2s(cpu_s.conf.current().vel(i))
                  << " gpu=" << math::v2s(gpu_s.conf.current().vel(i)) << std::endl;
        ++errors;
      }
      const math::Vec cf_diff = cpu_s.conf.old().constraint_force(i) -
                                 gpu_s.conf.old().constraint_force(i);
      if (math::abs(cf_diff) > cf_atol + cf_rtol * math::abs(cpu_s.conf.old().constraint_force(i))) {
        std::cerr << label << ": constraint_force mismatch at atom " << i
                  << ": cpu=" << math::v2s(cpu_s.conf.old().constraint_force(i))
                  << " gpu=" << math::v2s(gpu_s.conf.old().constraint_force(i)) << std::endl;
        ++errors;
      }
    }

    for (unsigned a = 0; a < 3; ++a) {
      for (unsigned b = 0; b < 3; ++b) {
        const double cpu_v = cpu_s.conf.old().virial_tensor(a, b);
        const double gpu_v = gpu_s.conf.old().virial_tensor(a, b);
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

  return run_case(stopo, sconf, sinput, "settle_gpu", quiet);
}
