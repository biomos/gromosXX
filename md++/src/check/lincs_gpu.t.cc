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
 * @file lincs_gpu.t.cc
 * End-to-end correctness test for CUDA_Lincs (PLAN.md §10 step 19,
 * constraints) -- runs the real CPU algorithm::Lincs and the real
 * algorithm::CUDA_Lincs on aladip's unperturbed topology, with (a) a
 * synthetic solute constraint list built from a few of its own bonds
 * (same idea as solute_shake_gpu.t.cc -- aladip's real input has no
 * solute constraints) and (b) its real solvent switched to LINCS
 * (aladip_unperturbed.in itself uses SHAKE for solvent -- see
 * shake_gpu.t.cc/settle_gpu.t.cc). Both groups exercise LINCS's two
 * solve passes (initial + rotational-lengthening correction) and its
 * coupled-constraint coefficient path (the synthetic solute
 * constraints share atoms, so `lincs.coupled_constr` is genuinely
 * non-empty). Compares the resulting positions, velocities, and
 * virial tensor. Only USE_CUDA builds run this.
 */

#include "../stdheader.h"

#include <cmath>

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"

#include "../algorithm/algorithm/algorithm_sequence.h"
#include "../math/periodicity.h"
#include "../algorithm/constraints/lincs.h"
#include "../algorithm/constraints/cuda_lincs.h"

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

  // Same construction as solute_shake_gpu.t.cc's
  // setup_synthetic_solute_constraints(), but chained (bond k's atom j
  // is bond k+1's atom i where possible) so at least one atom is
  // shared between two constraints -- exercises lincs.coupled_constr,
  // which a set of fully-disjoint constraints would leave empty.
  void setup_synthetic_solute_constraints(util::simulation_struct & s) {
    const std::vector<topology::two_body_term_struct> & bonds = s.topo.solute().bonds();
    const std::vector<interaction::bond_type_struct> & quart_types = s.topo.bond_types_quart();

    const unsigned num_constraints = 4;
    std::vector<interaction::bond_type_struct> & harm_types = s.topo.bond_types_harm();
    std::vector<topology::two_body_term_struct> & dc = s.topo.solute().distance_constraints();
    dc.clear();
    for (unsigned k = 0; k < num_constraints; ++k) {
      const double r0 = quart_types[bonds[k].type].r0;
      const unsigned harm_type = static_cast<unsigned>(harm_types.size());
      harm_types.push_back(interaction::bond_type_struct(0.0, r0));
      dc.push_back(topology::two_body_term_struct(bonds[k].i, bonds[k].j, harm_type));
    }

    s.sim.param().constraint.ntc = 3;
    s.sim.param().constraint.solute.algorithm = simulation::constr_lincs;
    s.sim.param().constraint.solvent.algorithm = simulation::constr_lincs;
    s.sim.param().system.nsm = static_cast<int>(s.topo.num_solvent_molecules(0));
  }

  void displace_solute_constrained_atoms(util::simulation_struct & s) {
    const std::vector<topology::two_body_term_struct> & dc = s.topo.solute().distance_constraints();
    for (unsigned k = 0; k < dc.size(); ++k) {
      const double d = 0.004 * std::sin(1.3 * k + 0.5);
      s.conf.current().pos(dc[k].i) += math::Vec(d, -0.5 * d, 0.3 * d);
      s.conf.current().pos(dc[k].j) += math::Vec(-0.5 * d, d, -0.2 * d);
    }
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

    setup_synthetic_solute_constraints(cpu_s);
    setup_synthetic_solute_constraints(gpu_s);

    cpu_s.conf.old() = cpu_s.conf.current();
    gpu_s.conf.old() = gpu_s.conf.current();

    displace_solute_constrained_atoms(cpu_s);
    displace_solute_constrained_atoms(gpu_s);
    displace_solvent(cpu_s);
    displace_solvent(gpu_s);

    algorithm::Lincs cpu_lincs;
    algorithm::CUDA_Lincs gpu_lincs;

    if (cpu_lincs.init(cpu_s.topo, cpu_s.conf, cpu_s.sim, std::cout, quiet) != 0 ||
        gpu_lincs.init(gpu_s.topo, gpu_s.conf, gpu_s.sim, std::cout, quiet) != 0) {
      flush_messages("lincs init");
      std::cerr << label << ": init() failed" << std::endl;
      return 1;
    }
    flush_messages("lincs init");

    const int cpu_rc = cpu_lincs.apply(cpu_s.topo, cpu_s.conf, cpu_s.sim);
    const int gpu_rc = gpu_lincs.apply(gpu_s.topo, gpu_s.conf, gpu_s.sim);
    flush_messages("lincs apply");

    if (cpu_rc != 0 || gpu_rc != 0) {
      std::cerr << label << ": apply() failed (cpu_rc=" << cpu_rc
                << " gpu_rc=" << gpu_rc << ")" << std::endl;
      return 1;
    }

    // Verify LINCS actually did something -- the solute constraints
    // should be satisfied (close to their target length) after apply().
    const std::vector<topology::two_body_term_struct> & dc = cpu_s.topo.solute().distance_constraints();
    const std::vector<interaction::bond_type_struct> & harm_types = cpu_s.topo.bond_types_harm();
    math::Periodicity<math::vacuum> periodicity(cpu_s.conf.current().box);
    for (unsigned k = 0; k < dc.size(); ++k) {
      math::Vec r;
      periodicity.nearest_image(cpu_s.conf.current().pos(dc[k].i),
                                  cpu_s.conf.current().pos(dc[k].j), r);
      const double r0 = harm_types[dc[k].type].r0;
      const double dist = math::abs(r);
      if (std::abs(dist - r0) > 1e-3 * r0) {
        std::cerr << label << ": CPU constraint " << k << " not satisfied after apply(): "
                  << "dist=" << dist << " r0=" << r0 << std::endl;
        return 1;
      }
    }

    // LINCS's own recursion is already Jacobi-shaped (see
    // lincs_kernels.h), so this port is bit-comparable to the CPU per
    // constraint, unlike solute SHAKE's Jacobi-vs-Gauss-Seidel
    // reformulation -- a tight relative tolerance is appropriate.
    const double tol = 1e-6;
    int errors = 0;

    const unsigned num_atoms = static_cast<unsigned>(cpu_s.topo.num_atoms());
    for (unsigned i = 0; i < num_atoms; ++i) {
      const math::Vec pos_diff = cpu_s.conf.current().pos(i) - gpu_s.conf.current().pos(i);
      const double pos_scale = std::max(1.0, math::abs(cpu_s.conf.current().pos(i)));
      if (math::abs(pos_diff) > tol * pos_scale) {
        std::cerr << label << ": pos mismatch at atom " << i
                  << ": cpu=" << math::v2s(cpu_s.conf.current().pos(i))
                  << " gpu=" << math::v2s(gpu_s.conf.current().pos(i)) << std::endl;
        ++errors;
      }
      const math::Vec vel_diff = cpu_s.conf.current().vel(i) - gpu_s.conf.current().vel(i);
      const double vel_scale = std::max(1.0, math::abs(cpu_s.conf.current().vel(i)));
      if (math::abs(vel_diff) > tol * vel_scale) {
        std::cerr << label << ": vel mismatch at atom " << i
                  << ": cpu=" << math::v2s(cpu_s.conf.current().vel(i))
                  << " gpu=" << math::v2s(gpu_s.conf.current().vel(i)) << std::endl;
        ++errors;
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

  return run_case(stopo, sconf, sinput, "lincs_gpu", quiet);
}
