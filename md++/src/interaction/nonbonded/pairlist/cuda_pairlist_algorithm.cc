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

  m_impl.init(topo, conf, sim, os, quiet);

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

  // prepare_cog() itself branches on atomic_cutoff (TILE_PAIRLIST_DESIGN.md
  // §4.2/step 7): always refreshes the GPU position mirror, and only does
  // the chargegroup cog/cell box-wrap when chargegroup-cutoff mode is
  // active. Safe to call unconditionally, and safe to call every step even
  // if atomic_cutoff is toggled between calls -- reorder()/build_candidates()/
  // classify_tiles() below re-read sim.param().pairlist.atomic_cutoff fresh
  // on every call too, so there's no stale chargegroup-mode state to fall
  // out of sync with.
  m_impl.prepare_cog(conf, topo, sim);

  return 0;
}

/**
 * Runs on the pairlist.skip_step cadence (see Nonbonded_Set::calculate_
 * interactions' pairlist_update check) -- the classification work, per
 * TILE_PAIRLIST_DESIGN.md §3 steps 1-5. m_impl.tiles() now holds real,
 * exclusion-checked, short/long-classified atom-pair tiles after this
 * call, consumed for real by compute_forces_energies() (TILE_PAIRLIST_
 * DESIGN.md §8/§9), called separately by CUDA_Nonbonded_Interaction --
 * not through this method's own `pairlist` parameter, which stays an
 * unused, empty interaction::PairlistContainer: the base
 * Pairlist_Algorithm::update() interface is CPU-Pairlist-shaped and
 * nothing reads this parameter for the CUDA path (kept only because
 * update() is virtual and must be implemented with this signature).
 *
 * The (expensive) candidate rebuild (reorder()+build_candidates()) is
 * decoupled from classification as of TILE_PAIRLIST_DESIGN.md §10 part
 * 2: m_impl.needs_candidate_rebuild() decides, via the skin-buffer
 * Verlet criterion, whether a rebuild is *also* due this
 * classification cycle, or whether the existing (possibly
 * several-classification-cycles-old, but still within the skin buffer)
 * candidate set is still safe to classify against. At skin == 0.0 this
 * always rebuilds, exactly as before this decoupling existed.
 */
void interaction::CUDA_Pairlist_Algorithm::update(topology::Topology & topo,
                                      configuration::Configuration & conf,
                                      simulation::Simulation &sim,
                                      interaction::PairlistContainer &pairlist,
                                      unsigned int begin, unsigned int end,
                                      unsigned int stride) {
  DEBUG(0, "cuda pairlist algorithm : update");

  // reorder()/build_candidates()/classify_tiles() each branch on
  // sim.param().pairlist.atomic_cutoff themselves (TILE_PAIRLIST_DESIGN.md
  // §4.2/step 7) -- chargegroup-cutoff and atomic-cutoff are both fully
  // supported now, so there's nothing to guard here.
  if (m_impl.needs_candidate_rebuild(conf, sim)) {
    m_impl.rebuild_candidates(conf, topo, sim);
  }
  m_impl.classify_tiles(conf, topo, sim);

  pairlist.clear();
}
