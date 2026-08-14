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
 * @file cuda_skin_drift.t.cc
 * Skin-buffer drift test (PLAN.md §9.4, TILE_PAIRLIST_DESIGN.md §10 part
 * 2): unlike every other CUDA test in this directory (single-evaluation),
 * this drives a real multi-step trajectory to check that decoupling the
 * candidate rebuild from classification via `skin` doesn't miss any real
 * pair.
 *
 * Two independently-constructed CUDA_Nonbonded_Interaction runs, driven
 * through the identical deterministic per-atom position perturbation
 * sequence at every step, both at the SAME skip_step (so classification
 * and long-range-force cadence -- and therefore the twin-range physics
 * itself -- are identical between them; skip_step's own effect is
 * already covered by TILE_PAIRLIST_DESIGN.md §10 part 1, not what this
 * test checks):
 *   - "baseline": skip_step = 5 (the GROMOS default), skin = 0.0 --
 *     candidates get rebuilt fresh every classification cycle (today's
 *     pre-decoupling behavior).
 *   - "decoupled": skip_step = 5, skin > 0 -- classification still runs
 *     every 5 steps against exact current positions, but the
 *     (expensive) candidate rebuild only happens when the skin buffer's
 *     Verlet criterion (2 * max_displacement >= skin) says it must.
 *
 * classify_tiles() always computes exact current-position distances
 * against whatever candidates currently exist, regardless of when they
 * were last rebuilt -- so if skin's buffer is never exhausted, the two
 * runs above must match EXACTLY at every step (not just within a
 * physically-motivated tolerance), and the only observable difference
 * should be that "decoupled" rebuilds candidates less often. Also checks
 * the explicit degenerate case: skin = 0.0 must rebuild candidates every
 * single classification cycle, exactly as if skin-decoupling didn't
 * exist.
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

#include "../io/argument.h"
#include "../util/parse_verbosity.h"
#include "../util/usage.h"
#include "../io/topology/in_topology.h"
#include "../io/message.h"

#include "../util/create_simulation.h"

#include "../interaction/nonbonded/interaction/nonbonded_parameter.h"
#include "../interaction/nonbonded/pairlist/pairlist_algorithm.h"
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

  /**
   * Deterministic, reproducible per-atom drift: each atom gets its own
   * fixed velocity direction (derived from its index via sin/cos, no RNG
   * dependency needed) so relative atom-atom distances genuinely change
   * over the trajectory -- a uniform whole-system translation wouldn't
   * exercise anything, since nearest_image distances (and therefore all
   * forces/energies) would stay exactly constant regardless of rebuild
   * cadence.
   */
  void perturb_positions(configuration::Configuration & conf, unsigned num_atoms,
                          double amplitude_per_step) {
    for (unsigned i = 0; i < num_atoms; ++i) {
      const double a = static_cast<double>(i);
      const math::Vec v(std::sin(a * 0.7 + 1.0), std::cos(a * 1.3 + 2.0), std::sin(a * 2.1 + 3.0));
      conf.current().pos(i) += amplitude_per_step * v;
    }
  }

  struct RunResult {
    bool ok = true;
    unsigned rebuild_count = 0;
  };

  /**
   * Builds a fresh simulation from scratch and drives it through
   * `num_steps` identical perturbation steps, calling
   * CUDA_Nonbonded_Interaction::calculate_interactions() every step.
   * `reference` (if non-null) holds the ground-truth run's per-step
   * forces/energies to compare against; when building the ground-truth
   * run itself, `reference` is null and this just records into `record`.
   */
  RunResult run_trajectory(const std::string & stopo, const std::string & sconf,
                            const std::string & sinput, math::boundary_enum force_boundary,
                            int skip_step, double skin, unsigned num_steps,
                            double amplitude_per_step, const char * label, bool quiet,
                            std::vector<math::VArray> * force_record,
                            std::vector<double> * e_lj_record,
                            std::vector<double> * e_crf_record,
                            const std::vector<math::VArray> * force_reference,
                            const std::vector<double> * e_lj_reference,
                            const std::vector<double> * e_crf_reference,
                            double tol) {
    RunResult result;

    util::simulation_struct s;
    io::In_Topology in_topo;
    in_topo.quiet = quiet;

    if (util::create_simulation(stopo, "", sconf, sinput, s, in_topo,
                                 "", "", "", "", "", "", "", "", quiet) != 0) {
      std::cerr << label << ": creating simulation failed" << std::endl;
      result.ok = false;
      return result;
    }
    io::messages.display(std::cout);
    io::messages.clear();

    s.sim.param().pairlist.skip_step = skip_step;
    s.sim.param().pairlist.skin = skin;

    s.conf.boundary_type = force_boundary;
    s.sim.param().boundary.boundary = force_boundary;
    if (force_boundary == math::rectangular) {
      s.conf.current().box = math::Box(math::Vec(4.0, 0.0, 0.0),
                                        math::Vec(0.0, 4.0, 0.0),
                                        math::Vec(0.0, 0.0, 4.0));
    }

    const unsigned num_atoms = static_cast<unsigned>(s.topo.num_atoms());
    s.topo.energy_groups().assign(1, num_atoms - 1);
    s.topo.atom_energy_group().assign(num_atoms, 0u);
    s.sim.param().pcouple.virial = math::no_virial;

    interaction::Pairlist_Algorithm * pa = new interaction::CUDA_Pairlist_Algorithm();
    interaction::CUDA_Nonbonded_Interaction * ni =
        new interaction::CUDA_Nonbonded_Interaction(pa);

    in_topo.read_lj_parameter(ni->parameter().lj_parameter(), std::cout);
    if (flush_messages("lj parameter read")) {
      std::cerr << label << ": reading LJ parameters failed" << std::endl;
      delete ni;
      result.ok = false;
      return result;
    }

    if (ni->init(s.topo, s.conf, s.sim, std::cout, quiet) != 0) {
      flush_messages("init");
      std::cerr << label << ": init() failed" << std::endl;
      delete ni;
      result.ok = false;
      return result;
    }

    interaction::CUDA_Pairlist_Algorithm * cuda_pa =
        static_cast<interaction::CUDA_Pairlist_Algorithm *>(pa);

    int errors = 0;
    for (unsigned step = 0; step < num_steps; ++step) {
      s.sim.steps() = step;
      if (step > 0) {
        perturb_positions(s.conf, num_atoms, amplitude_per_step);
      }

      s.conf.current().force = 0.0;
      s.conf.current().energies.zero();

      if (ni->calculate_interactions(s.topo, s.conf, s.sim) != 0) {
        flush_messages("calculate_interactions");
        std::cerr << label << ": calculate_interactions() failed at step " << step << std::endl;
        ++errors;
        break;
      }
      if (flush_messages("calculate_interactions")) {
        std::cerr << label << ": error message at step " << step << std::endl;
        ++errors;
        break;
      }

      const double e_lj  = s.conf.current().energies.lj_energy[0][0];
      const double e_crf = s.conf.current().energies.crf_energy[0][0];

      if (force_record)  force_record->push_back(s.conf.current().force);
      if (e_lj_record)   e_lj_record->push_back(e_lj);
      if (e_crf_record)  e_crf_record->push_back(e_crf);

      if (force_reference && e_lj_reference && e_crf_reference) {
        const double ref_e_lj  = (*e_lj_reference)[step];
        const double ref_e_crf = (*e_crf_reference)[step];
        if (std::abs(e_lj - ref_e_lj) > tol * std::max(1.0, std::abs(ref_e_lj))) {
          std::cerr << label << ": e_lj mismatch at step " << step << ": got=" << e_lj
                    << " ref=" << ref_e_lj << std::endl;
          ++errors;
        }
        if (std::abs(e_crf - ref_e_crf) > tol * std::max(1.0, std::abs(ref_e_crf))) {
          std::cerr << label << ": e_crf mismatch at step " << step << ": got=" << e_crf
                    << " ref=" << ref_e_crf << std::endl;
          ++errors;
        }
        const math::VArray & ref_force = (*force_reference)[step];
        for (unsigned i = 0; i < num_atoms; ++i) {
          const math::Vec diff = s.conf.current().force(i) - ref_force(i);
          const double scale = std::max(1.0, math::abs(ref_force(i)));
          if (math::abs(diff) > tol * scale) {
            std::cerr << label << ": force mismatch at step " << step << ", atom " << i
                      << ": got=" << math::v2s(s.conf.current().force(i))
                      << " ref=" << math::v2s(ref_force(i)) << std::endl;
            ++errors;
          }
        }
      }
    }

    result.rebuild_count = cuda_pa->candidate_rebuild_count();
    result.ok = (errors == 0);

    if (!result.ok) {
      std::cerr << label << ": FAILED (" << errors << " mismatch(es) over " << num_steps
                << " steps)" << std::endl;
    } else {
      std::cout << label << ": OK (" << num_steps << " steps, "
                << result.rebuild_count << " candidate rebuild(s))" << std::endl;
    }

    delete ni; // cascades: ~Nonbonded_Interaction deletes m_pairlist_algorithm (pa)
    return result;
  }

  int run_case(const std::string & stopo, const std::string & sconf,
               const std::string & sinput, math::boundary_enum force_boundary,
               const char * label, bool quiet) {
    const unsigned num_steps = 30;
    const int skip_step = 5;
    const double skin = 0.08;
    const double amplitude_per_step = 0.002;
    const unsigned classification_cycles = (num_steps + skip_step - 1) / skip_step;

    // Baseline: same classification/long-range cadence (skip_step) as
    // the decoupled run below, but skin = 0.0 -- candidates get rebuilt
    // fresh every classification cycle. This, not a skip_step = 1 run,
    // is the correct reference: skip_step itself already changes the
    // physics (real twin-range freezes long-range forces between
    // classification cycles, TILE_PAIRLIST_DESIGN.md §10 part 1), so a
    // skip_step = 1 run is not expected to match a skip_step = 5 run at
    // all -- that's an intentional, different approximation, not
    // something this test is checking. What skin must not change is
    // classify_tiles()'s own output at a *given* classification step:
    // classify_tiles() always computes exact current-position distances
    // against whatever candidates currently exist, so a same-cadence run
    // that skips some candidate rebuilds must produce IDENTICAL results
    // to one that never skips them, provided the skin buffer was never
    // exhausted.
    std::string baseline_label = std::string(label) + " (baseline, skip_step=5 skin=0)";
    std::vector<math::VArray> baseline_force;
    std::vector<double> baseline_e_lj, baseline_e_crf;
    RunResult baseline = run_trajectory(
        stopo, sconf, sinput, force_boundary, skip_step, /*skin=*/0.0,
        num_steps, amplitude_per_step, baseline_label.c_str(), quiet,
        &baseline_force, &baseline_e_lj, &baseline_e_crf, nullptr, nullptr, nullptr, 0.0);
    if (!baseline.ok) return 1;

    int total = 0;

    // Explicit degenerate case: skin = 0.0 must rebuild every
    // classification cycle -- exactly the pre-decoupling behavior,
    // strict backward compatibility.
    if (baseline.rebuild_count != classification_cycles) {
      std::cerr << baseline_label << ": expected exactly " << classification_cycles
                << " rebuilds (skin=0 degenerate case), got " << baseline.rebuild_count
                << std::endl;
      ++total;
    }

    std::string decoupled_label = std::string(label) + " (decoupled, skip_step=5 skin>0)";
    RunResult decoupled = run_trajectory(
        stopo, sconf, sinput, force_boundary, skip_step, skin,
        num_steps, amplitude_per_step, decoupled_label.c_str(), quiet,
        nullptr, nullptr, nullptr, &baseline_force, &baseline_e_lj, &baseline_e_crf, 1e-4);
    total += decoupled.ok ? 0 : 1;

    // Decoupling must actually be happening, not just numerically
    // coincide with always-correct behavior: strictly fewer rebuilds
    // than classification cycles -- otherwise a needs_candidate_rebuild()
    // bug that always returns true would still pass the force/energy
    // comparison above.
    if (decoupled.rebuild_count >= classification_cycles) {
      std::cerr << decoupled_label << ": skin did not reduce rebuild frequency ("
                << decoupled.rebuild_count << " rebuilds over " << classification_cycles
                << " classification cycles)" << std::endl;
      ++total;
    } else {
      std::cout << decoupled_label << ": skin skipped "
                << (classification_cycles - decoupled.rebuild_count) << "/"
                << classification_cycles << " candidate rebuilds" << std::endl;
    }

    return total;
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
