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
 * @file special_force_after_gpu_bonded.t.cc
 * Regression test for a real bug: a CPU-only "special" force term
 * (position restraints, distance restraints, etc. -- every class under
 * interaction/special/) that runs *after* GPU-native bonded/nonbonded
 * terms inside the same Forcefield must see the GPU-computed force
 * already published to conf.current().force before it does its own
 * `+=`, and that addition must survive the eventual GPU-mirror publish
 * (which *overwrites* conf.current().force from the mirror, not
 * merges). Before Interaction::needs_fresh_cpu_force() defaulted to
 * true (interaction.h), a CPU-only special force running after
 * CUDA_Quartic_Bond_Interaction would silently add onto a stale/zero
 * array, and that addition would later be wiped out when the mirror's
 * force was published -- the bonded contribution simply vanished.
 * Reported by a real user run (POSITIONRES active alongside GPU-native
 * bonded/nonbonded), surfacing as escalating LINCS "too much rotation"
 * warnings over thousands of steps -- exactly what dropped force from
 * a chunk of the system looks like. This test reproduces the mechanism
 * directly (a synthetic special-like force term, not real position
 * restraints, so it needs no @posresspec/@refpos input files) rather
 * than the full multi-thousand-step scenario.
 */

#include "../stdheader.h"

#include <cmath>

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"

#include "../algorithm/algorithm/algorithm_sequence.h"
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

  // Stand-in for any interaction/special/* class: adds a fixed,
  // per-atom force directly onto conf.current().force -- the exact
  // shape of Position_Restraint_Interaction/Distance_Restraint_
  // Interaction/etc (a plain CPU array reference, `+=`), without
  // needing real restraint-specification input files.
  class Synthetic_Special_Force : public interaction::Interaction {
  public:
    Synthetic_Special_Force() : Interaction("SyntheticSpecialForce") {}

    int init(topology::Topology &, configuration::Configuration &,
             simulation::Simulation &, std::ostream &, bool) override {
      return 0;
    }

    int calculate_interactions(topology::Topology & topo,
                                configuration::Configuration & conf,
                                simulation::Simulation &) override {
      const unsigned num_atoms = static_cast<unsigned>(topo.num_atoms());
      for (unsigned i = 0; i < num_atoms; ++i)
        conf.current().force(i) += math::Vec(1.0, 2.0, 3.0);
      return 0;
    }
  };

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

    // Forcefield::~Forcefield() deletes every pointer it holds.
    interaction::Forcefield cpu_ff, gpu_ff;
    cpu_ff.push_back(new interaction::Quartic_Bond_Interaction());
    cpu_ff.push_back(new Synthetic_Special_Force());
    gpu_ff.push_back(new interaction::CUDA_Quartic_Bond_Interaction());
    gpu_ff.push_back(new Synthetic_Special_Force());

    if (cpu_ff.init(cpu_s.topo, cpu_s.conf, cpu_s.sim, std::cout, quiet) != 0 ||
        gpu_ff.init(gpu_s.topo, gpu_s.conf, gpu_s.sim, std::cout, quiet) != 0) {
      flush_messages("forcefield init");
      std::cerr << label << ": init() failed" << std::endl;
      return 1;
    }
    flush_messages("forcefield init");

    // calculate_interactions() itself zeroes force (Forcefield's own
    // convention) -- no manual zero needed here.
    if (cpu_ff.calculate_interactions(cpu_s.topo, cpu_s.conf, cpu_s.sim) != 0 ||
        gpu_ff.calculate_interactions(gpu_s.topo, gpu_s.conf, gpu_s.sim) != 0) {
      flush_messages("calculate_interactions");
      std::cerr << label << ": calculate_interactions() failed" << std::endl;
      return 1;
    }
    flush_messages("calculate_interactions");

    const double tol = 5e-3; // see quartic_bond_gpu.t.cc's tolerance comment
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

    if (errors) {
      std::cerr << label << ": FAILED (" << errors << " mismatch(es))" << std::endl;
    } else {
      std::cout << label << ": OK (bonded + synthetic special force both present, "
                << num_atoms << " atoms)" << std::endl;
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

  return run_case(stopo, sconf, sinput, "special_force_after_gpu_bonded", quiet);
}
