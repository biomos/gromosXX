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
 * @file shake_gpu.t.cc
 * End-to-end correctness test for CUDA_Shake (PLAN.md §10 step 15,
 * constraints) -- runs the real CPU algorithm::Shake and the real
 * algorithm::CUDA_Shake on aladip's unperturbed topology/configuration
 * (NTC=1: solvent-only constraints, exactly CUDA_Shake's v1 scope) and
 * compares the resulting positions, velocities, constraint forces, and
 * virial tensor.
 *
 * aladip_unperturbed.in's own starting configuration already satisfies
 * every solvent constraint almost exactly (SHAKE would converge in 0-1
 * iterations, barely exercising the algorithm), so this test displaces
 * conf.current().pos for solvent atoms by a small deterministic offset
 * before shaking -- conf.old().pos (the reference geometry SHAKE
 * corrects back towards) is left untouched, exactly the shape of a real
 * MD step's post-integration, pre-constraint state. Only USE_CUDA
 * builds run this.
 */

#include "../stdheader.h"

#include <cmath>

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"

#include "../algorithm/algorithm/algorithm_sequence.h"
#include "../math/periodicity.h"
#include "../algorithm/constraints/shake.h"
#include "../algorithm/constraints/cuda_shake.h"
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

  // Small deterministic per-component displacement, distinct per atom
  // so different constraints within a molecule are perturbed
  // differently -- large enough to force several SHAKE iterations
  // (aladip's solvent is SPC water, O-H/H-H bond lengths ~0.1-0.16 nm),
  // small enough to stay well inside SHAKE's basin of convergence.
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

    if (cpu_s.sim.param().constraint.ntc != 1) {
      std::cerr << label << ": expected NTC == 1 (solvent-only constraints) "
                << "for this topology/input" << std::endl;
      return 1;
    }

    // util::create_simulation() doesn't populate conf.old() from the
    // configuration file (only conf.current() has real data) -- set
    // old() = current() first, matching a real MD step's starting
    // point, then displace only current() (the "free", unconstrained
    // positions SHAKE corrects back towards old()).
    cpu_s.conf.old() = cpu_s.conf.current();
    gpu_s.conf.old() = gpu_s.conf.current();

    displace_solvent(cpu_s);
    displace_solvent(gpu_s);

    algorithm::Shake cpu_shake(cpu_s.sim.param().constraint.solute.shake_tolerance,
                               cpu_s.sim.param().constraint.solvent.shake_tolerance);
    algorithm::CUDA_Shake gpu_shake(gpu_s.sim.param().constraint.solute.shake_tolerance,
                                    gpu_s.sim.param().constraint.solvent.shake_tolerance);

    if (cpu_shake.init(cpu_s.topo, cpu_s.conf, cpu_s.sim, std::cout, quiet) != 0 ||
        gpu_shake.init(gpu_s.topo, gpu_s.conf, gpu_s.sim, std::cout, quiet) != 0) {
      flush_messages("shake init");
      std::cerr << label << ": init() failed" << std::endl;
      return 1;
    }
    flush_messages("shake init");

    const int cpu_rc = cpu_shake.apply(cpu_s.topo, cpu_s.conf, cpu_s.sim);
    const int gpu_rc = gpu_shake.apply(gpu_s.topo, gpu_s.conf, gpu_s.sim);
    flush_messages("shake apply");

    // CUDA_Shake leaves pos/vel GPU-resident (mark_gpu_dirty(), no eager
    // CPU publish) -- a real consumer only pays for the publish when it
    // actually needs the value (Algorithm_Sequence::run()'s generic
    // flush, or program/md.cc's output-cadence flush); this test reads
    // gpu_s.conf directly, so it must request that publish explicitly,
    // same as any other direct consumer would.
    gpu_s.sim.cuda().flush_gpu_dirty(gpu_s.conf, gpu::MIRROR_POS | gpu::MIRROR_VEL | gpu::MIRROR_CONSTRAINT_FORCE);

    if (cpu_rc != 0 || gpu_rc != 0) {
      std::cerr << label << ": apply() failed (cpu_rc=" << cpu_rc
                << " gpu_rc=" << gpu_rc << ")" << std::endl;
      return 1;
    }

    // Position stays tight -- GPU now computes in FPL_TYPE (float under
    // FP_PRECISION 1/2), but position itself still agrees with CPU
    // (double) to within this tolerance. constraint_force/vel are
    // position differences divided by dt^2/dt, which amplifies ordinary
    // float32-vs-double rounding noise by orders of magnitude (~1/dt^2
    // for constraint_force) -- verified this is exactly that (not a
    // bug) by checking pos itself agrees tightly. See m_shake_gpu.t.cc's
    // fuller comment.
    const double tol = 1e-6;
    const double cf_atol = 1.0, cf_rtol = 2e-3;
    const double vel_atol = 2e-3, vel_rtol = 5e-4;
    int errors = 0;

    const unsigned num_atoms = static_cast<unsigned>(cpu_s.topo.num_atoms());
    const unsigned first_solvent = static_cast<unsigned>(cpu_s.topo.num_solute_atoms());
    for (unsigned i = first_solvent; i < num_atoms; ++i) {
      const math::Vec pos_diff = cpu_s.conf.current().pos(i) - gpu_s.conf.current().pos(i);
      if (math::abs(pos_diff) > tol) {
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

  return run_case(stopo, sconf, sinput, "shake_gpu", quiet);
}
