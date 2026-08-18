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
 * @file cuda_nonbonded_interaction.t.cc
 * End-to-end correctness test for CUDA_Nonbonded_Interaction (PLAN.md
 * §10 step 9, TILE_PAIRLIST_DESIGN.md §8/§9) -- unlike
 * pairlist_cuda_equivalence.t.cc (pairlist only) and
 * lj_crf_tile_kernel.t.cc (kernel only, hand-built tile), this drives
 * the real CUDA_Pairlist_Algorithm + CUDA_Nonbonded_Interaction pair
 * exactly as create_nonbonded.cc wires them for accelerator = cuda, and
 * compares the resulting per-atom forces and per-energy-group-pair
 * LJ/CRF energy matrix against a direct CPU reference built from
 * Standard_Pairlist_Algorithm's real pairlist (all four buckets:
 * solute/solvent x short/long) summed through
 * interaction::Nonbonded_Term::lj_crf_interaction -- the same primitive
 * the CPU innerloop uses.
 *
 * aladip's own test input defines 2 energy groups (NEGR = 2), used as-is
 * here (no override) -- both sides read topo.atom_energy_group() for
 * their energy-group-pair bucketing (CUDA_Nonbonded_Interaction's tile
 * kernel now buckets by [eg_i][eg_j] instead of a single flat total), so
 * this is a real multi-energy-group correctness test, not a
 * single-group-forced one. Only USE_CUDA builds run this.
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

#include "../math/periodicity.h"
#include "../interaction/nonbonded/interaction/nonbonded_parameter.h"
#include "../interaction/nonbonded/interaction/nonbonded_term.h"
#include "../interaction/nonbonded/pairlist/pairlist.h"
#include "../interaction/nonbonded/pairlist/pairlist_algorithm.h"
#include "../interaction/nonbonded/pairlist/standard_pairlist_algorithm.h"
#include "../interaction/nonbonded/pairlist/cuda_pairlist_algorithm.h"
#include "../interaction/nonbonded/interaction/cuda_nonbonded_interaction.h"

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

  struct Reference {
    std::vector<math::Vec> force;
    // Per-energy-group-pair matrices, matching
    // configuration::Energy::lj_energy/crf_energy's [gi][gj] shape.
    std::vector<std::vector<double> > lj_energy;
    std::vector<std::vector<double> > crf_energy;
  };

  /**
   * Real CPU pairlist (Standard_Pairlist_Algorithm, all four buckets) +
   * direct Nonbonded_Term::lj_crf_interaction summation, bucketed by
   * energy-group pair exactly like nonbonded_innerloop.cc's CPU inner
   * loop (storage.energies.lj_energy[eg_i][eg_j] += e_lj) -- the CPU-side
   * ground truth this test compares CUDA_Nonbonded_Interaction against.
   */
  template <math::boundary_enum B>
  Reference compute_reference(topology::Topology & topo,
                               configuration::Configuration & conf,
                               simulation::Simulation & sim,
                               interaction::Nonbonded_Parameter & param) {
    interaction::Standard_Pairlist_Algorithm cpu_pa;
    cpu_pa.init(topo, conf, sim, std::cout, true);
    cpu_pa.prepare(topo, conf, sim);
    interaction::PairlistContainer pl;
    pl.resize(static_cast<unsigned>(topo.num_atoms()));
    cpu_pa.update(topo, conf, sim, pl, 0, static_cast<unsigned>(topo.num_atoms()), 1);

    interaction::Nonbonded_Term term;
    term.init(sim);
    math::Periodicity<B> periodicity(conf.current().box);

    const unsigned num_groups = static_cast<unsigned>(topo.energy_groups().size());

    Reference ref;
    ref.force.resize(topo.num_atoms(), math::Vec(0.0, 0.0, 0.0));
    ref.lj_energy.assign(num_groups, std::vector<double>(num_groups, 0.0));
    ref.crf_energy.assign(num_groups, std::vector<double>(num_groups, 0.0));

    auto accumulate = [&](interaction::Pairlist & bucket) {
      for (unsigned i = 0; i < bucket.size(); ++i) {
        for (unsigned int j : bucket[i]) {
          math::Vec r;
          periodicity.nearest_image(conf.current().pos(i), conf.current().pos(j), r);
          const interaction::lj_parameter_struct & lj =
              param.lj_parameter(topo.iac(i), topo.iac(j));
          const double q = topo.charge(i) * topo.charge(j);
          double f = 0.0, e_lj = 0.0, e_crf = 0.0;
          term.lj_crf_interaction(r, lj.c6, lj.c12, q, f, e_lj, e_crf);
          ref.force[i] += f * r;
          ref.force[j] -= f * r;
          const unsigned eg_i = topo.atom_energy_group(i);
          const unsigned eg_j = topo.atom_energy_group(j);
          ref.lj_energy[eg_i][eg_j]  += e_lj;
          ref.crf_energy[eg_i][eg_j] += e_crf;
        }
      }
    };
    accumulate(pl.solute_short);
    accumulate(pl.solute_long);
    accumulate(pl.solvent_short);
    accumulate(pl.solvent_long);

    // 1,4 pairs ("LJ exceptions"): excluded from the regular pairlist
    // above (all_exclusion = exclusion UNION one_four_pair), added back
    // in with their own scaled cs6/cs12 LJ parameters and Coulomb-scaled
    // CRF -- matches Nonbonded_Outerloop::one_four_outerloop's
    // unconditional call (nonbonded_set.cc) exactly, now that
    // gpu::launch_one_four (cuda_pairlist_algorithm_impl.cu) implements
    // this on the GPU side too.
    for (unsigned i = 0; i < topo.num_solute_atoms(); ++i) {
      for (unsigned int j : topo.one_four_pair(i)) {
        math::Vec r;
        periodicity.nearest_image(conf.current().pos(i), conf.current().pos(j), r);
        const interaction::lj_parameter_struct & lj =
            param.lj_parameter(topo.iac(i), topo.iac(j));
        const double q = topo.charge(i) * topo.charge(j);
        double f = 0.0, e_lj = 0.0, e_crf = 0.0;
        term.lj_crf_interaction(r, lj.cs6, lj.cs12, q, f, e_lj, e_crf, 0,
                                 param.get_coulomb_scaling());
        ref.force[i] += f * r;
        ref.force[j] -= f * r;
        const unsigned eg_i = topo.atom_energy_group(i);
        const unsigned eg_j = topo.atom_energy_group(j);
        ref.lj_energy[eg_i][eg_j]  += e_lj;
        ref.crf_energy[eg_i][eg_j] += e_crf;
      }
    }

    return ref;
  }

  int run_case(const std::string & stopo, const std::string & sconf,
               const std::string & sinput, math::boundary_enum force_boundary,
               const char * label, bool quiet) {
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

    // Exact match required, same as pairlist_cuda_equivalence.t.cc.
    s.sim.param().pairlist.skin = 0.0;

    // rf_excluded defaults to true (parameter.h's "new standard") and is
    // now genuinely implemented on GPU (gpu::launch_rf_excluded) -- this
    // test's CPU reference is deliberately just the regular pairlist-
    // driven LJ/CRF sum, not that separate contribution, so disable it
    // here to keep comparing like for like. rf_excluded_gpu.t.cc is the
    // dedicated test for that term.
    s.sim.param().nonbonded.rf_excluded = false;

    s.conf.boundary_type = force_boundary;
    s.sim.param().boundary.boundary = force_boundary;
    if (force_boundary == math::rectangular) {
      s.conf.current().box = math::Box(math::Vec(4.0, 0.0, 0.0),
                                        math::Vec(0.0, 4.0, 0.0),
                                        math::Vec(0.0, 0.0, 4.0));
    }

    const unsigned num_atoms = static_cast<unsigned>(s.topo.num_atoms());

    // v1 scope: no virial (CUDA_Nonbonded_Interaction::init() gate).
    s.sim.param().pcouple.virial = math::no_virial;

    interaction::Nonbonded_Parameter param;
    in_topo.read_lj_parameter(param.lj_parameter(), std::cout);
    if (flush_messages("lj parameter read")) {
      std::cerr << label << ": reading LJ parameters failed" << std::endl;
      return 1;
    }

    Reference ref;
    if (force_boundary == math::vacuum) {
      ref = compute_reference<math::vacuum>(s.topo, s.conf, s.sim, param);
    } else {
      ref = compute_reference<math::rectangular>(s.topo, s.conf, s.sim, param);
    }

    // GPU side: exactly the pairing create_nonbonded.cc uses for
    // accelerator = cuda.
    interaction::Pairlist_Algorithm * pa = new interaction::CUDA_Pairlist_Algorithm();
    interaction::CUDA_Nonbonded_Interaction * ni =
        new interaction::CUDA_Nonbonded_Interaction(pa);

    in_topo.read_lj_parameter(ni->parameter().lj_parameter(), std::cout);
    if (flush_messages("lj parameter read (gpu)")) {
      std::cerr << label << ": reading LJ parameters failed (gpu)" << std::endl;
      delete ni;
      return 1;
    }

    if (ni->init(s.topo, s.conf, s.sim, std::cout, quiet) != 0) {
      flush_messages("gpu init");
      std::cerr << label << ": CUDA_Nonbonded_Interaction::init() failed" << std::endl;
      delete ni;
      return 1;
    }

    s.conf.current().force = 0.0;
    s.conf.current().energies.zero();

    if (ni->calculate_interactions(s.topo, s.conf, s.sim) != 0) {
      flush_messages("gpu calculate_interactions");
      std::cerr << label << ": CUDA_Nonbonded_Interaction::calculate_interactions() failed"
                << std::endl;
      delete ni;
      return 1;
    }
    if (flush_messages("gpu calculate_interactions")) {
      delete ni;
      return 1;
    }

    const unsigned num_groups = static_cast<unsigned>(s.topo.energy_groups().size());

    int errors = 0;
    // FPL_TYPE is float in the default (mixed-precision) build.
    const double tol = 1e-4;

    for (unsigned gi = 0; gi < num_groups; ++gi) {
      for (unsigned gj = 0; gj < num_groups; ++gj) {
        const double gpu_e_lj  = s.conf.current().energies.lj_energy[gi][gj];
        const double gpu_e_crf = s.conf.current().energies.crf_energy[gi][gj];
        const double ref_e_lj  = ref.lj_energy[gi][gj];
        const double ref_e_crf = ref.crf_energy[gi][gj];
        if (std::abs(gpu_e_lj - ref_e_lj) > tol * std::max(1.0, std::abs(ref_e_lj))) {
          std::cerr << label << ": e_lj mismatch at group (" << gi << "," << gj
                    << "): gpu=" << gpu_e_lj << " cpu=" << ref_e_lj << std::endl;
          ++errors;
        }
        if (std::abs(gpu_e_crf - ref_e_crf) > tol * std::max(1.0, std::abs(ref_e_crf))) {
          std::cerr << label << ": e_crf mismatch at group (" << gi << "," << gj
                    << "): gpu=" << gpu_e_crf << " cpu=" << ref_e_crf << std::endl;
          ++errors;
        }
      }
    }
    for (unsigned i = 0; i < num_atoms; ++i) {
      const math::Vec diff = s.conf.current().force(i) - ref.force[i];
      const double scale = std::max(1.0, math::abs(ref.force[i]));
      if (math::abs(diff) > tol * scale) {
        std::cerr << label << ": force mismatch at atom " << i
                  << ": gpu=" << math::v2s(s.conf.current().force(i))
                  << " cpu=" << math::v2s(ref.force[i]) << std::endl;
        ++errors;
      }
    }

    if (errors) {
      std::cerr << label << ": FAILED (" << errors << " mismatch(es))" << std::endl;
    } else {
      std::cout << label << ": OK (" << num_atoms << " atoms, " << num_groups
                << " energy group(s))" << std::endl;
    }

    delete ni; // cascades: ~Nonbonded_Interaction deletes m_pairlist_algorithm (pa)
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
  total += run_case(stopo, sconf, sinput, math::vacuum, "vacuum", quiet);
  total += run_case(stopo, sconf, sinput, math::rectangular, "rectangular", quiet);

  return total;
}
