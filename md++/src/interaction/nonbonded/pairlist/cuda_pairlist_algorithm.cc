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
 * @file cuda_pairlist_algorithm.cc
 * CUDA accelerated pairlist algorithm
 */

#include "stdheader.h"

#include "algorithm/algorithm.h"
#include "topology/topology.h"
#include "simulation/simulation.h"
#include "configuration/configuration.h"

#include "math/periodicity.h"

#include "pairlist.h"
#include "pairlist_algorithm.h"
#include "cuda_pairlist_algorithm.h"

#include "util/debug.h"
#include "util/template_split.h"

#include "gpu/cuda/utils.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist

interaction::CUDA_Pairlist_Algorithm::CUDA_Pairlist_Algorithm()
          : Pairlist_Algorithm()
{}

int interaction::CUDA_Pairlist_Algorithm::init(
                                      topology::Topology & topo,
                                      configuration::Configuration & conf,
                                      simulation::Simulation & sim,
                                      std::ostream & os,
                                      bool quiet)
{
  const math::boundary_enum b = conf.boundary_type;
  if (b != math::vacuum && b != math::rectangular) {
    io::messages.add(
      "CUDA_Pairlist_Algorithm only supports vacuum and rectangular "
      "boundary conditions so far (TILE_PAIRLIST_DESIGN.md §4.1); "
      "triclinic/truncoct cell lists are real future work, not yet done.",
      "CUDA_Pairlist_Algorithm", io::message::error);
    return 1;
  }

  // TILE_PAIRLIST_DESIGN.md §4.2: chargegroup-cutoff first, atomic-cutoff
  // as a separate follow-up step. Calling with atomic_cutoff=true today
  // would silently skip the chargegroup cog/cell build this pairlist
  // currently depends on entirely, producing an empty (wrong) candidate
  // list rather than a working atomic-cutoff one -- hard-error instead.
  if (sim.param().pairlist.atomic_cutoff) {
    io::messages.add(
      "CUDA_Pairlist_Algorithm does not support atomic_cutoff yet "
      "(TILE_PAIRLIST_DESIGN.md §4.2: chargegroup-cutoff lands first).",
      "CUDA_Pairlist_Algorithm", io::message::error);
    return 1;
  }

  if (!quiet)
    os << "\tcuda pairlist algorithm\n";
  return 0;
}

/**
 * Runs every step (see Nonbonded_Interaction::calculate_interactions --
 * "shared memory do this only once [per step]"), unlike update() below,
 * which only runs on the pairlist.skip_step cadence. Box-wrapping and cog
 * computation must happen every step since atoms move every step; the
 * expensive block/candidate build must not (TILE_PAIRLIST_DESIGN.md §3
 * step 1: "expensive, O(N), infrequent").
 */
int interaction::CUDA_Pairlist_Algorithm::prepare(
                                      topology::Topology & topo,
                                      configuration::Configuration & conf,
                                      simulation::Simulation & sim)
{
  DEBUG(0, "cuda pairlist algorithm : prepare");

  m_impl.set_cutoff(sim.param().pairlist.cutoff_short,
	     sim.param().pairlist.cutoff_long);

  // init()'s atomic_cutoff hard-error is a one-time, construction-time
  // check; sim.param().pairlist.atomic_cutoff can still be flipped at
  // runtime after init() ran (confirmed by check_forcefield.cc's "atomic
  // cutoff" comparison test, which reuses one already-initialized
  // Forcefield/CUDA_Pairlist_Algorithm under both settings) -- prepare_cog
  // already no-ops when atomic_cutoff is true, so this stays a no-op too,
  // but don't rely on init() alone to keep that guarantee. See update()'s
  // matching guard for why this matters beyond just prepare_cog: reorder()/
  // build_candidates() must not run against a chargegroup-mode cog/cell
  // build that was silently skipped this call.
  if (!sim.param().pairlist.atomic_cutoff) {
    m_impl.prepare_cog(conf, topo, sim);
  }

  return 0;
}

/**
 * Runs on the pairlist.skip_step cadence (see Nonbonded_Set::calculate_
 * interactions' pairlist_update check) -- the real candidate-build work,
 * per TILE_PAIRLIST_DESIGN.md §3 steps 2-4. Still ends by leaving the
 * CPU-facing `pairlist` (interaction::PairlistContainer) an explicit,
 * warned-about dummy: nothing downstream (no force kernel, no
 * classification pass) consumes m_impl.tiles() yet, so there is nothing
 * real to report through this container -- see PAIRLIST_PLAN.md §5(A).
 */
void interaction::CUDA_Pairlist_Algorithm::update(topology::Topology & topo,
                                      configuration::Configuration & conf,
                                      simulation::Simulation &sim,
                                      interaction::PairlistContainer &pairlist,
                                      unsigned int begin, unsigned int end,
                                      unsigned int stride) {
  DEBUG(0, "cuda pairlist algorithm : update");

  // Same runtime-toggle concern as prepare()'s guard above: only run
  // reorder()/build_candidates() when this call's atomic_cutoff setting
  // actually matches what prepare_cog() built this cycle. Observed for
  // real: check_forcefield.cc's "atomic cutoff" check flips
  // sim.param().pairlist.atomic_cutoff to true on an already-initialized
  // CUDA_Pairlist_Algorithm and calls update() again -- without this
  // guard, reorder() ran thrust::sort_by_key while the CUDA context was
  // apparently left in a bad state by that surrounding test sequence
  // (manifested as "invalid device ordinal" from CUB/Thrust) instead of
  // cleanly no-op'ing like every other part of this class already does
  // for atomic_cutoff.
  if (!sim.param().pairlist.atomic_cutoff) {
    m_impl.reorder(conf, topo, sim);
    m_impl.build_candidates(conf, topo, sim);
  } else if (!m_warned_atomic_cutoff) {
    io::messages.add(
      "CUDA_Pairlist_Algorithm does not support atomic_cutoff yet "
      "(TILE_PAIRLIST_DESIGN.md §4.2); skipping the candidate build this "
      "call.",
      "CUDA_Pairlist_Algorithm", io::message::warning);
    m_warned_atomic_cutoff = true;
  }

  // TODO(cleanup): dummy placeholder, see PAIRLIST_PLAN.md §5(A). Not the
  // real force-consumable output -- produces an intentionally empty
  // pairlist (zero nonbonded pairs) so that selecting accelerator=cuda
  // fails loudly (visibly wrong, zero energy, plus this warning) rather
  // than quietly running a plausible-looking but fake result.
  pairlist.clear();
  if (!m_warned_dummy) {
    io::messages.add(
      "CUDA pairlist algorithm is a placeholder (PAIRLIST_PLAN.md); "
      "it produces zero nonbonded pairs until the real tile-based GPU "
      "pairlist (TILE_PAIRLIST_DESIGN.md) is implemented.",
      "CUDA_Pairlist_Algorithm", io::message::warning);
    m_warned_dummy = true;
  }
}
