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
 * @file pressure_barostat_gpu.t.cc
 * Pressure_Calculation (plain host-side, O(9) matrix code, no CUDA
 * kernel of its own -- narrows gpu_mirror_touches() to MIRROR_VIRIAL,
 * see pressure_calculation.h) and Berendsen_Barostat<gpuBackend> (O(9)
 * mu/box computation on the host, O(num_atoms) position scaling as a
 * real GPU kernel, GPU-resident -- see berendsen_barostat_gpu.cc) are
 * both exercised here through a real Algorithm_Sequence::run(), and
 * compared against the identical CPU-only sequence:
 * CUDA_Quartic_Bond_Interaction (GPU force calc) -> Pressure_Calculation
 * -> Berendsen_Barostat<gpuBackend> -> CUDA_Quartic_Bond_Interaction
 * again. The second force calc must see the barostat's GPU-scaled
 * positions (whether published to CPU yet or not -- configuration_view()
 * resolves that transparently), not a stale pre-barostat copy. Only
 * USE_CUDA builds run this.
 */

#include "../stdheader.h"

#include <cmath>

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"

#include "../algorithm/algorithm/algorithm_sequence.h"
#include "../algorithm/pressure/pressure_calculation.h"
#include "../algorithm/pressure/berendsen_barostat.h"
#include "../interaction/interaction.h"
#include "../interaction/forcefield/forcefield.h"
#include "../interaction/bonded/quartic_bond_interaction.h"
#include "../interaction/bonded/cuda_quartic_bond_interaction.h"

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

    // Enable isotropic pressure coupling so Berendsen_Barostat actually
    // scales the box/positions (aladip_unperturbed.in itself has no
    // pressure coupling block worth relying on).
    cpu_s.sim.param().pcouple.calculate = true;
    cpu_s.sim.param().pcouple.scale = math::pcouple_isotropic;
    cpu_s.sim.param().pcouple.virial = math::atomic_virial;
    cpu_s.sim.param().pcouple.compressibility = 4.575e-4;
    cpu_s.sim.param().pcouple.tau = 0.5;
    cpu_s.sim.param().pcouple.pres0 = math::Matrix(0.06102);
    gpu_s.sim.param() = cpu_s.sim.param();

    // Forcefield::~Forcefield() deletes every pointer it holds -- these
    // must be heap-allocated, not stack locals.
    interaction::Forcefield cpu_ff, gpu_ff;
    cpu_ff.push_back(new interaction::Quartic_Bond_Interaction());
    gpu_ff.push_back(new interaction::CUDA_Quartic_Bond_Interaction());

    algorithm::Pressure_Calculation cpu_pcalc, gpu_pcalc;
    algorithm::Berendsen_Barostat<util::cpuBackend> cpu_baro;
    algorithm::Berendsen_Barostat<util::gpuBackend> gpu_baro;

    if (cpu_ff.init(cpu_s.topo, cpu_s.conf, cpu_s.sim, std::cout, quiet) != 0 ||
        gpu_ff.init(gpu_s.topo, gpu_s.conf, gpu_s.sim, std::cout, quiet) != 0) {
      flush_messages("forcefield init");
      std::cerr << label << ": forcefield init() failed" << std::endl;
      return 1;
    }
    flush_messages("forcefield init");

    // clean=false: these sequences hold pointers to stack locals
    // (cpu_ff/gpu_ff/*pcalc/*baro below) -- Algorithm_Sequence::
    // ~Algorithm_Sequence() deletes its contents unless told not to.
    algorithm::Algorithm_Sequence cpu_force_seq(false), gpu_force_seq(false);
    cpu_force_seq.push_back(&cpu_ff);
    gpu_force_seq.push_back(&gpu_ff);

    algorithm::Algorithm_Sequence cpu_pressure_seq(false), gpu_pressure_seq(false);
    cpu_pressure_seq.push_back(&cpu_pcalc);
    cpu_pressure_seq.push_back(&cpu_baro);
    gpu_pressure_seq.push_back(&gpu_pcalc);
    gpu_pressure_seq.push_back(&gpu_baro);

    // 1st GPU force calc -- populates conf.current().virial_tensor.
    if (cpu_force_seq.run(cpu_s.topo, cpu_s.conf, cpu_s.sim) != 0 ||
        gpu_force_seq.run(gpu_s.topo, gpu_s.conf, gpu_s.sim) != 0) {
      flush_messages("first force calc");
      std::cerr << label << ": first force calculation failed" << std::endl;
      return 1;
    }
    flush_messages("first force calc");

    // gpu_ff's bonded/nonbonded terms atomicAdd their virial straight
    // into the GPU mirror's virial_tensor (no CPU round trip) --
    // gpu_force_seq.run() only contains gpu_ff itself, so nothing in
    // that sequence triggers the publish; must flush explicitly before
    // reading gpu_s.conf.current().virial_tensor below.
    gpu_s.sim.cuda().flush_gpu_dirty(gpu_s.conf, gpu::MIRROR_VIRIAL);

    // Pressure_Calculation reads conf.old(), not conf.current() -- feed
    // it the virial this step just computed plus a fixed synthetic
    // kinetic energy tensor (identical on both sides; nothing in this
    // minimal standalone sequence otherwise populates it).
    cpu_s.conf.old().virial_tensor = cpu_s.conf.current().virial_tensor;
    gpu_s.conf.old().virial_tensor = gpu_s.conf.current().virial_tensor;
    cpu_s.conf.old().kinetic_energy_tensor = math::Matrix(100.0);
    gpu_s.conf.old().kinetic_energy_tensor = math::Matrix(100.0);

    // Pressure_Calculation + Berendsen_Barostat -- CPU-only, no CUDA
    // awareness; Algorithm_Sequence::run()'s default gpu_mirror_touches()
    // handling (PLAN.md §10 step 16) is what's actually under test here.
    if (cpu_pressure_seq.run(cpu_s.topo, cpu_s.conf, cpu_s.sim) != 0 ||
        gpu_pressure_seq.run(gpu_s.topo, gpu_s.conf, gpu_s.sim) != 0) {
      flush_messages("pressure/barostat");
      std::cerr << label << ": pressure/barostat failed" << std::endl;
      return 1;
    }
    flush_messages("pressure/barostat");

    // 2nd GPU force calc -- must see the barostat's host-scaled
    // positions, not a stale pre-barostat GPU-resident copy.
    if (cpu_force_seq.run(cpu_s.topo, cpu_s.conf, cpu_s.sim) != 0 ||
        gpu_force_seq.run(gpu_s.topo, gpu_s.conf, gpu_s.sim) != 0) {
      flush_messages("second force calc");
      std::cerr << label << ": second force calculation failed" << std::endl;
      return 1;
    }
    flush_messages("second force calc");

    // gpu_ff's CUDA_Quartic_Bond_Interaction writes force directly into
    // the GPU mirror and mark_gpu_dirty()s it -- gpu_force_seq.run()
    // publishes it before the *next* algorithm in that same sequence
    // (there isn't one here), so this standalone comparison needs its
    // own explicit publish (see angle_gpu.t.cc's comment). Same for
    // Berendsen_Barostat<gpuBackend>'s POS write, deferred all the way
    // through the second force calc (nothing in gpu_pressure_seq/
    // gpu_force_seq ever needed a CPU-fresh copy) -- see leap_frog_gpu.
    // t.cc's identical comment.
    gpu_s.sim.cuda().flush_gpu_dirty(gpu_s.conf, gpu::MIRROR_FORCE | gpu::MIRROR_POS);

    const double tol = 5e-3; // see quartic_bond_gpu.t.cc's tolerance comment
    // Berendsen_Barostat<gpuBackend>'s position-scaling kernel runs in
    // FPL_TYPE (float under FP_PRECISION=2) -- mu itself is computed on
    // the host in double and matches the CPU path exactly (see the box
    // comparison below, still 1e-9), but pos(i) = mu * pos(i) loses
    // precision to ~1e-7 in the multiply/store. 1e-9 (appropriate when
    // this class was CPU-only on both sides) is too tight for that.
    const double pos_tol = 1e-5;
    int errors = 0;

    const unsigned num_atoms = static_cast<unsigned>(cpu_s.topo.num_atoms());
    for (unsigned i = 0; i < num_atoms; ++i) {
      const math::Vec pos_diff = cpu_s.conf.current().pos(i) - gpu_s.conf.current().pos(i);
      const double pos_scale = std::max(1.0, math::abs(cpu_s.conf.current().pos(i)));
      if (math::abs(pos_diff) > pos_tol * pos_scale) {
        std::cerr << label << ": post-barostat pos mismatch at atom " << i
                  << ": cpu=" << math::v2s(cpu_s.conf.current().pos(i))
                  << " gpu=" << math::v2s(gpu_s.conf.current().pos(i)) << std::endl;
        ++errors;
      }
      const math::Vec force_diff = cpu_s.conf.current().force(i) - gpu_s.conf.current().force(i);
      const double scale = std::max(1.0, math::abs(cpu_s.conf.current().force(i)));
      if (math::abs(force_diff) > tol * scale) {
        std::cerr << label << ": post-barostat 2nd force-calc mismatch at atom " << i
                  << ": cpu=" << math::v2s(cpu_s.conf.current().force(i))
                  << " gpu=" << math::v2s(gpu_s.conf.current().force(i)) << std::endl;
        ++errors;
      }
    }

    for (unsigned a = 0; a < 3; ++a) {
      for (unsigned b = 0; b < 3; ++b) {
        const double cpu_box = cpu_s.conf.current().box(a)(b);
        const double gpu_box = gpu_s.conf.current().box(a)(b);
        if (std::abs(cpu_box - gpu_box) > 1e-9) {
          std::cerr << label << ": box(" << a << "," << b << ") mismatch: cpu="
                    << cpu_box << " gpu=" << gpu_box << std::endl;
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

  return run_case(stopo, sconf, sinput, "pressure_barostat_gpu", quiet);
}
