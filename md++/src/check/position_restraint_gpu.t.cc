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
 * @file position_restraint_gpu.t.cc
 * End-to-end correctness test for CUDA_Position_Restraint_Interaction --
 * runs the real CPU Position_Restraint_Interaction and the real
 * CUDA_Position_Restraint_Interaction on aladip's unperturbed topology/
 * configuration and compares the resulting per-atom forces and
 * per-energy-group posrest energies. aladip's own input has no restraint
 * specification, so this test builds a small synthetic
 * topo.position_restraints() + conf.special().reference_positions/
 * bfactors by hand (same idea cuda_nonbonded_interaction.t.cc uses for
 * synthetic energy-group overrides), covering both posrest_on and
 * posrest_bfactor modes. Only USE_CUDA builds run this.
 */

#include "../stdheader.h"

#include <cmath>

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"

#include "../algorithm/algorithm/algorithm_sequence.h"
#include "../interaction/interaction.h"
#include "../interaction/special/position_restraint_interaction.h"
#include "../interaction/special/cuda_position_restraint_interaction.h"

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

  // Restrain every solute atom to a fixed offset from its current
  // position (guarantees a nonzero restraint force for every term,
  // without needing a real @refpos/@posresspec file) -- identical setup
  // applied to both the cpu_s and gpu_s simulation_structs.
  void build_synthetic_restraints(util::simulation_struct & s,
                                   bool use_bfactor) {
    const unsigned num_atoms = static_cast<unsigned>(s.topo.num_atoms());

    s.conf.special().reference_positions.resize(num_atoms);
    s.conf.special().bfactors.resize(num_atoms);

    s.topo.position_restraints().clear();
    for (unsigned i = 0; i < num_atoms; ++i) {
      s.topo.position_restraints().push_back(topology::position_restraint_struct(i));
      s.conf.special().reference_positions(i) =
          s.conf.current().pos(i) + math::Vec(0.01, -0.02, 0.03);
      // Varying, always-positive bfactor -- only consumed under
      // posrest_bfactor, harmless (unread) otherwise.
      s.conf.special().bfactors(i) = 0.5 + 0.1 * (i % 5);
    }

    s.sim.param().posrest.posrest =
        use_bfactor ? simulation::posrest_bfactor : simulation::posrest_on;
    s.sim.param().posrest.force_constant = 2.5e4;
  }

  int run_case(const std::string & stopo, const std::string & sconf,
               const std::string & sinput, bool use_bfactor,
               const char * label, bool quiet) {
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

    build_synthetic_restraints(cpu_s, use_bfactor);
    build_synthetic_restraints(gpu_s, use_bfactor);

    // Both sides need a zeroed force/energy starting point, same as
    // Forcefield::calculate_interactions() provides in a real run --
    // both Interaction classes accumulate (+=), don't overwrite.
    cpu_s.conf.current().force = 0.0;
    gpu_s.conf.current().force = 0.0;
    for (auto & e : cpu_s.conf.current().energies.posrest_energy) e = 0.0;
    for (auto & e : gpu_s.conf.current().energies.posrest_energy) e = 0.0;

    interaction::Position_Restraint_Interaction cpu_pr;
    interaction::CUDA_Position_Restraint_Interaction gpu_pr;

    if (cpu_pr.init(cpu_s.topo, cpu_s.conf, cpu_s.sim, std::cout, quiet) != 0 ||
        gpu_pr.init(gpu_s.topo, gpu_s.conf, gpu_s.sim, std::cout, quiet) != 0) {
      flush_messages("position_restraint init");
      std::cerr << label << ": init() failed" << std::endl;
      return 1;
    }
    flush_messages("position_restraint init");

    if (cpu_pr.calculate_interactions(cpu_s.topo, cpu_s.conf, cpu_s.sim) != 0 ||
        gpu_pr.calculate_interactions(gpu_s.topo, gpu_s.conf, gpu_s.sim) != 0) {
      flush_messages("position_restraint calculate_interactions");
      std::cerr << label << ": calculate_interactions() failed" << std::endl;
      return 1;
    }
    flush_messages("position_restraint calculate_interactions");

    // CUDA_Position_Restraint_Interaction writes force directly into the
    // GPU mirror and mark_gpu_dirty()s it (see angle_gpu.t.cc's comment
    // for why this standalone test needs its own explicit publish).
    gpu_s.sim.cuda().flush_gpu_dirty(gpu_s.conf, gpu::MIRROR_FORCE);
    gpu_s.sim.cuda().flush_gpu_dirty(gpu_s.conf, gpu::MIRROR_ENERGY);

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

    const unsigned num_groups =
        static_cast<unsigned>(cpu_s.conf.current().energies.posrest_energy.size());
    for (unsigned g = 0; g < num_groups; ++g) {
      const double cpu_e = cpu_s.conf.current().energies.posrest_energy[g];
      const double gpu_e = gpu_s.conf.old().energies.posrest_energy[g];
      const double escale = std::max(1.0, std::abs(cpu_e));
      if (std::abs(cpu_e - gpu_e) > tol * escale) {
        std::cerr << label << ": posrest energy mismatch in group " << g
                  << ": cpu=" << cpu_e << " gpu=" << gpu_e << std::endl;
        ++errors;
      }
    }

    if (errors) {
      std::cerr << label << ": FAILED (" << errors << " mismatch(es))" << std::endl;
    } else {
      std::cout << label << ": OK (" << num_atoms << " restrained atoms, "
                << (use_bfactor ? "posrest_bfactor" : "posrest_on") << ")" << std::endl;
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

  int result = 0;
  result += run_case(stopo, sconf, sinput, false, "position_restraint_gpu (posrest_on)", quiet);
  result += run_case(stopo, sconf, sinput, true, "position_restraint_gpu (posrest_bfactor)", quiet);
  return result;
}
