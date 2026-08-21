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
 * @file lattice_shift_gpu.t.cc
 * Validates Lattice_Shift_Tracker<gpuBackend> against the CPU reference
 * (math::Periodicity<b>::put_chargegroups_into_box_saving_shifts(),
 * math/periodicity.cc). aladip's own configuration already has every
 * chargegroup inside the box (nothing to wrap, a trivially-passing but
 * uninteresting test), so this deliberately displaces every atom by a
 * whole box vector first -- guaranteeing every chargegroup crosses the
 * boundary and must be wrapped, exercising both the position rewrap and
 * the lattice_shifts bookkeeping this GPU port is actually for.
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
#include "../algorithm/integration/lattice_shift.h"

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

    const unsigned num_atoms = static_cast<unsigned>(cpu_s.topo.num_atoms());

    // Displace every atom by exactly one box vector (0) -- guarantees
    // every chargegroup's centre of geometry/first-atom position is
    // now outside the box on that axis and must be wrapped back in,
    // with a non-zero shift recorded.
    for (util::simulation_struct * s : {&cpu_s, &gpu_s}) {
      const math::Vec offset = s->conf.current().box(0);
      for (unsigned i = 0; i < num_atoms; ++i)
        s->conf.current().pos(i) += offset;
      s->conf.special().lattice_shifts = 0.0;
    }

    algorithm::Lattice_Shift_Tracker<util::cpuBackend> cpu_lst;
    algorithm::Lattice_Shift_Tracker<util::gpuBackend> gpu_lst;

    if (cpu_lst.init(cpu_s.topo, cpu_s.conf, cpu_s.sim, std::cout, quiet) != 0 ||
        gpu_lst.init(gpu_s.topo, gpu_s.conf, gpu_s.sim, std::cout, quiet) != 0) {
      flush_messages("init");
      std::cerr << label << ": init() failed" << std::endl;
      return 1;
    }
    flush_messages("init");

    if (cpu_lst.apply(cpu_s.topo, cpu_s.conf, cpu_s.sim) != 0 ||
        gpu_lst.apply(gpu_s.topo, gpu_s.conf, gpu_s.sim) != 0) {
      flush_messages("apply");
      std::cerr << label << ": apply() failed" << std::endl;
      return 1;
    }
    flush_messages("apply");

    // apply() defers the publish (mark_gpu_dirty(), no CPU round trip)
    // -- this test reads conf.current().pos/conf.special().
    // lattice_shifts directly, so it must ask for the publish itself,
    // same as leap_frog_gpu.t.cc's identical fix.
    gpu_s.sim.cuda().flush_gpu_dirty(gpu_s.conf, gpu::MIRROR_POS | gpu::MIRROR_LATTICE_SHIFT);

    const double tol = 1e-4;
    int errors = 0;

    for (unsigned i = 0; i < num_atoms; ++i) {
      const math::Vec dx = cpu_s.conf.current().pos(i) - gpu_s.conf.current().pos(i);
      const double xscale = std::max(1.0, math::abs(cpu_s.conf.current().pos(i)));
      if (math::abs(dx) > tol * xscale) {
        std::cerr << label << ": pos mismatch at atom " << i
                  << ": cpu=" << math::v2s(cpu_s.conf.current().pos(i))
                  << " gpu=" << math::v2s(gpu_s.conf.current().pos(i)) << std::endl;
        ++errors;
      }
      const math::Vec ds = cpu_s.conf.special().lattice_shifts(i) - gpu_s.conf.special().lattice_shifts(i);
      const double sscale = std::max(1.0, math::abs(cpu_s.conf.special().lattice_shifts(i)));
      if (math::abs(ds) > tol * sscale) {
        std::cerr << label << ": lattice_shift mismatch at atom " << i
                  << ": cpu=" << math::v2s(cpu_s.conf.special().lattice_shifts(i))
                  << " gpu=" << math::v2s(gpu_s.conf.special().lattice_shifts(i)) << std::endl;
        ++errors;
      }
    }

    // Sanity check the test itself actually exercised wrapping (not a
    // vacuously-passing all-zero comparison): at least one atom's
    // shift must be non-zero, given the whole-box displacement above.
    bool any_nonzero = false;
    for (unsigned i = 0; i < num_atoms; ++i) {
      if (math::abs(cpu_s.conf.special().lattice_shifts(i)) > tol) { any_nonzero = true; break; }
    }
    if (!any_nonzero) {
      std::cerr << label << ": test setup did not actually trigger any wrapping "
                << "(all lattice_shifts are zero) -- CPU reference itself is suspect" << std::endl;
      ++errors;
    }

    if (errors == 0) {
      std::cout << label << ": OK (" << num_atoms
                << " atoms, CPU/GPU lattice-shift wrap+bookkeeping match)" << std::endl;
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

  return run_case(stopo, sconf, sinput, "lattice_shift_gpu", quiet);
}
