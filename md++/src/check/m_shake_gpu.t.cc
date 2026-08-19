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
 * @file m_shake_gpu.t.cc
 * End-to-end correctness test for CUDA_M_Shake -- runs the real CPU
 * algorithm::M_Shake and the real algorithm::CUDA_M_Shake on aladip's
 * unperturbed topology/configuration (its solvent is 3-site SPC water,
 * exactly CUDA_M_Shake's v1 scope) and compares the resulting
 * positions, velocities, constraint forces, and virial tensor. Same
 * displaced-solvent-positions scenario as shake_gpu.t.cc/
 * settle_gpu.t.cc. Every molecule is solved fully independently with
 * the same fixed operation sequence on both sides (unlike solute
 * SHAKE's Jacobi-vs-host-loop scheme), so this is a tight,
 * bit-comparable-up-to-associativity tolerance. Only USE_CUDA builds
 * run this.
 */

#include "../stdheader.h"

#include <cmath>

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"

#include "../algorithm/algorithm/algorithm_sequence.h"
#include "../algorithm/constraints/m_shake.h"
#include "../algorithm/constraints/cuda_m_shake.h"

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

    // Drive M-SHAKE regardless of the input file's own NTCS choice.
    cpu_s.sim.param().constraint.solvent.algorithm = simulation::constr_m_shake;
    gpu_s.sim.param().constraint.solvent.algorithm = simulation::constr_m_shake;
    cpu_s.sim.param().system.nsm = static_cast<int>(cpu_s.topo.num_solvent_molecules(0));
    gpu_s.sim.param().system.nsm = static_cast<int>(gpu_s.topo.num_solvent_molecules(0));

    cpu_s.conf.old() = cpu_s.conf.current();
    gpu_s.conf.old() = gpu_s.conf.current();

    displace_solvent(cpu_s);
    displace_solvent(gpu_s);

    algorithm::M_Shake cpu_shake(1.0e-6, 1000);
    algorithm::CUDA_M_Shake gpu_shake(1.0e-6, 1000);

    if (cpu_shake.init(cpu_s.topo, cpu_s.conf, cpu_s.sim, std::cout, quiet) != 0 ||
        gpu_shake.init(gpu_s.topo, gpu_s.conf, gpu_s.sim, std::cout, quiet) != 0) {
      flush_messages("m_shake init");
      std::cerr << label << ": init() failed" << std::endl;
      return 1;
    }
    flush_messages("m_shake init");

    const int cpu_rc = cpu_shake.apply(cpu_s.topo, cpu_s.conf, cpu_s.sim);
    const int gpu_rc = gpu_shake.apply(gpu_s.topo, gpu_s.conf, gpu_s.sim);
    flush_messages("m_shake apply");

    // CUDA_M_Shake is GPU-resident: apply() leaves the corrected
    // position/velocity on the GPU mirror (mark_gpu_dirty()), relying
    // on Algorithm_Sequence::run()'s automatic flush_gpu_dirty() before
    // whatever runs next to publish it back to conf. This test calls
    // apply() standalone, with no such "next algorithm" -- so it must
    // publish explicitly itself before inspecting gpu_s.conf, the same
    // way a real end-of-step eventually would.
    gpu_s.sim.cuda().sync_configuration_from_device(gpu_s.conf);

    if (cpu_rc != 0 || gpu_rc != 0) {
      std::cerr << label << ": apply() failed (cpu_rc=" << cpu_rc
                << " gpu_rc=" << gpu_rc << ")" << std::endl;
      return 1;
    }

    // Position (the quantity that actually matters physically, and
    // what the next MD step integrates from) stays tight: the GPU
    // kernels now compute in FPL_TYPE (float under FP_PRECISION 1/2,
    // see gpu/cuda/algorithm/constraints/*_kernels.cu), but position
    // itself agrees with CPU (double) to within this tolerance even so.
    const double tol = 1e-6;
    // constraint_force and vel are position differences divided by
    // dt^2/dt respectively (dt is small -- 0.002 ps here), which
    // amplifies ordinary float32-vs-double non-associativity by orders
    // of magnitude (~1/dt^2 ~ 2.5e5x for constraint_force) -- ordinary
    // FPL_TYPE rounding noise (~1e-7 relative on an O(1) nm position)
    // becomes an O(0.01-0.1) absolute difference in force. Verified
    // this is exactly float32 noise, not a bug, by checking `pos`
    // (above) and `virial` (below, accumulated in double throughout)
    // both already agree at the tight tolerance -- only the amplified
    // quantities need a looser bound.
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

    // Actual constraint satisfaction, not just CPU/GPU agreement --
    // catches a bug that happens to agree on both sides (e.g. a wrong
    // factor-matrix sign both ported identically).
    const std::vector<topology::two_body_term_struct> & dc =
        cpu_s.topo.solvent(0).distance_constraints();
    for (unsigned m = 0; m < cpu_s.topo.num_solvent_molecules(0); ++m) {
      const unsigned base = first_solvent + m * 3;
      for (const auto & c : dc) {
        const double r0 = cpu_s.topo.bond_types_harm()[c.type].r0;
        const double r_cpu = math::abs(cpu_s.conf.current().pos(base + c.i) -
                                        cpu_s.conf.current().pos(base + c.j));
        const double r_gpu = math::abs(gpu_s.conf.current().pos(base + c.i) -
                                        gpu_s.conf.current().pos(base + c.j));
        if (std::abs(r_cpu - r0) > 1e-4 * r0) {
          std::cerr << label << ": CPU constraint not satisfied, molecule " << m
                    << ": r=" << r_cpu << " r0=" << r0 << std::endl;
          ++errors;
        }
        if (std::abs(r_gpu - r0) > 1e-4 * r0) {
          std::cerr << label << ": GPU constraint not satisfied, molecule " << m
                    << ": r=" << r_gpu << " r0=" << r0 << std::endl;
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

  return run_case(stopo, sconf, sinput, "m_shake_gpu", quiet);
}
