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
 * @file pairlist_cuda_equivalence.t.cc
 * Pairlist equivalence test: CUDA_Pairlist_Algorithm vs
 * Standard_Pairlist_Algorithm. TILE_PAIRLIST_DESIGN.md §5/§6 -- the
 * actual correctness gate for the tile-based GPU pairlist work
 * (steps 3-5): asserts the two algorithms produce exactly the same
 * per-atom short/long pair sets (including exclusions) on the same
 * topology/configuration, at skin = 0. Only USE_CUDA builds run this
 * (CUDA_Pairlist_Algorithm doesn't exist otherwise).
 */

#include "../stdheader.h"

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"

#include "../algorithm/algorithm/algorithm_sequence.h"

#include "../io/argument.h"
#include "../util/parse_verbosity.h"
#include "../util/usage.h"
#include "../io/topology/in_topology.h"
#include "../io/message.h"

#include "../util/create_simulation.h"

#include "../interaction/nonbonded/pairlist/pairlist.h"
#include "../interaction/nonbonded/pairlist/pairlist_algorithm.h"
#include "../interaction/nonbonded/pairlist/standard_pairlist_algorithm.h"
#include "../interaction/nonbonded/pairlist/cuda_pairlist_algorithm.h"

#include "check.h"

#ifdef XXMPI
  #include <mpi.h>
#endif

namespace {

  // Displays and clears any io::messages queued since the last call,
  // returning true if any were error-or-worse severity. io::messages.add
  // only queues messages -- nothing shows them unless something calls
  // display() -- so skipping this would silently swallow real errors
  // (e.g. a candidate-capacity overflow) exactly like a passing run.
  bool flush_messages(const char * label) {
    const io::message::severity_enum sev = io::messages.display(std::cerr);
    io::messages.clear();
    return sev >= io::message::error;
  }

  interaction::PairlistContainer run_cpu(topology::Topology & topo,
                                          configuration::Configuration & conf,
                                          simulation::Simulation & sim) {
    interaction::Standard_Pairlist_Algorithm alg;
    alg.init(topo, conf, sim, std::cout, true);
    alg.prepare(topo, conf, sim);
    interaction::PairlistContainer pl;
    pl.resize(static_cast<unsigned>(topo.num_atoms()));
    alg.update(topo, conf, sim, pl, 0, static_cast<unsigned>(topo.num_atoms()), 1);
    flush_messages("cpu");
    return pl;
  }

  // Returns false (and leaves `out` untouched) if init() rejects the
  // configuration (e.g. an unsupported boundary type) or anything during
  // prepare()/update() reported an error (e.g. a candidate-capacity
  // overflow) -- the caller reports either as a hard failure, not a
  // silent skip.
  bool run_gpu(topology::Topology & topo,
               configuration::Configuration & conf,
               simulation::Simulation & sim,
               interaction::PairlistContainer & out) {
    interaction::CUDA_Pairlist_Algorithm alg;
    if (alg.init(topo, conf, sim, std::cout, true) != 0) {
      flush_messages("gpu init");
      return false;
    }
    alg.prepare(topo, conf, sim);
    // update()'s own output parameter is an intentional dummy
    // (PAIRLIST_PLAN.md §5(A)) -- the real result comes from the tiles
    // it builds internally, via to_pairlist_container() below.
    interaction::PairlistContainer dummy;
    dummy.resize(static_cast<unsigned>(topo.num_atoms()));
    alg.update(topo, conf, sim, dummy, 0, static_cast<unsigned>(topo.num_atoms()), 1);
    if (flush_messages("gpu update")) return false;
    out = alg.to_pairlist_container(topo);
    return true;
  }

  void sort_rows(interaction::PairlistContainer & pl) {
    interaction::Pairlist * buckets[4] = {
      &pl.solute_short, &pl.solute_long, &pl.solvent_short, &pl.solvent_long
    };
    for (interaction::Pairlist * bucket : buckets) {
      for (std::vector<unsigned int> & row : *bucket) {
        std::sort(row.begin(), row.end());
      }
    }
  }

  int compare_bucket(const char * name,
                      const interaction::Pairlist & cpu,
                      const interaction::Pairlist & gpu) {
    int errors = 0;
    const size_t n = std::min(cpu.size(), gpu.size());
    if (cpu.size() != gpu.size()) {
      std::cerr << "  " << name << ": size mismatch, cpu=" << cpu.size()
                << " gpu=" << gpu.size() << std::endl;
      ++errors;
    }
    for (size_t a = 0; a < n; ++a) {
      if (cpu[a] != gpu[a]) {
        std::cerr << "  " << name << " mismatch for atom " << a << ":\n"
                  << "    cpu:";
        for (unsigned j : cpu[a]) std::cerr << " " << j;
        std::cerr << "\n    gpu:";
        for (unsigned j : gpu[a]) std::cerr << " " << j;
        std::cerr << std::endl;
        ++errors;
      }
    }
    return errors;
  }

  int compare(const interaction::PairlistContainer & cpu,
              const interaction::PairlistContainer & gpu) {
    int errors = 0;
    errors += compare_bucket("solute_short",  cpu.solute_short,  gpu.solute_short);
    errors += compare_bucket("solute_long",   cpu.solute_long,   gpu.solute_long);
    errors += compare_bucket("solvent_short", cpu.solvent_short, gpu.solvent_short);
    errors += compare_bucket("solvent_long",  cpu.solvent_long,  gpu.solvent_long);
    return errors;
  }

  int run_case(const std::string & stopo, const std::string & sconf,
               const std::string & sinput, math::boundary_enum force_boundary,
               const char * label, bool quiet) {
    // Fresh simulation per case (not reused across cases): pairlist
    // prepare() mutates chargegroup positions in place (box-wrapping),
    // so reusing one topo/conf across differently-shaped boundary cases
    // would make each case depend on the previous one's side effects.
    util::simulation_struct s;
    io::In_Topology in_topo;
    in_topo.quiet = quiet;

    if (util::create_simulation(stopo, "", sconf, sinput, s, in_topo,
                                 "", "", "", "", "", "", "", "", quiet) != 0) {
      std::cerr << label << ": creating simulation failed" << std::endl;
      return 1;
    }
    io::messages.display(std::cout);
    io::messages.clear();

    // Exact match required (TILE_PAIRLIST_DESIGN.md §5): skin=0 means
    // the GPU candidate radius is exactly cutoff_long, same as the CPU
    // reference's rebuild radius.
    s.sim.param().pairlist.skin = 0.0;

    s.conf.boundary_type = force_boundary;
    s.sim.param().boundary.boundary = force_boundary;
    if (force_boundary == math::rectangular) {
      // Big enough to comfortably contain this (small) system with room
      // to spare -- the equivalence property doesn't depend on the exact
      // box size, only that both algorithms see the identical box.
      s.conf.current().box = math::Box(math::Vec(4.0, 0.0, 0.0),
                                        math::Vec(0.0, 4.0, 0.0),
                                        math::Vec(0.0, 0.0, 4.0));
    }

    interaction::PairlistContainer cpu_pl = run_cpu(s.topo, s.conf, s.sim);

    interaction::PairlistContainer gpu_pl;
    if (!run_gpu(s.topo, s.conf, s.sim, gpu_pl)) {
      std::cerr << label << ": CUDA_Pairlist_Algorithm rejected this "
                << "configuration or reported an error (see io::messages "
                << "output above)" << std::endl;
      return 1;
    }

    sort_rows(cpu_pl);
    sort_rows(gpu_pl);

    const int errors = compare(cpu_pl, gpu_pl);
    if (errors) {
      std::cerr << label << ": FAILED (" << errors << " mismatching row(s)/bucket(s))"
                << std::endl;
    } else {
      std::cout << label << ": OK ("
                << cpu_pl.solute_short.size() + cpu_pl.solvent_short.size()
                << " atoms compared)" << std::endl;
    }
    return errors;
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
  // aladip's own configuration is a truncated-octahedron box, which
  // neither the CPU reference used here for chargegroup-cutoff nor the
  // GPU pairlist attempt to handle identically for boundary-crossing
  // pairs -- boundary_type is overridden per case below precisely to
  // sidestep that and test vacuum/rectangular specifically (§4.1's v1
  // scope), reusing aladip's real solute+solvent topology/coordinates.
  total += run_case(stopo, sconf, sinput, math::vacuum, "vacuum", quiet);
  total += run_case(stopo, sconf, sinput, math::rectangular, "rectangular", quiet);

  return total;
}
