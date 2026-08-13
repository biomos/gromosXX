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

/**
 * calculate center of geometries
 */
int interaction::CUDA_Pairlist_Algorithm::prepare(
                                      topology::Topology & topo,
                                      configuration::Configuration & conf,
                                      simulation::Simulation & sim)
{
  DEBUG(0, "cuda pairlist algorithm : prepare");

  m_impl.set_cutoff(sim.param().pairlist.cutoff_short,
	     sim.param().pairlist.cutoff_long);

  if (!sim.param().pairlist.atomic_cutoff) {

    // first put the chargegroups into the box
    m_impl.prepare_cog(conf, topo, sim);
  } // chargegroup based cutoff

  // assign all cogs / atoms to grid cells and reorder
  m_impl.reorder(conf, topo, sim);

  return 0;

}

void interaction::CUDA_Pairlist_Algorithm::update(topology::Topology & topo,
                                      configuration::Configuration & conf,
                                      simulation::Simulation &sim,
                                      interaction::PairlistContainer &pairlist,
                                      unsigned int begin, unsigned int end,
                                      unsigned int stride) {
  DEBUG(0, "cuda pairlist algorithm : update");
  // TODO(cleanup): dummy placeholder, see PAIRLIST_PLAN.md §5(A)/§6 step 3.
  // Not the real tile-based GPU pairlist (TILE_PAIRLIST_DESIGN.md) --
  // produces an intentionally empty pairlist (zero nonbonded pairs) so
  // that selecting accelerator=cuda fails loudly (visibly wrong, zero
  // energy, plus this warning) rather than quietly running a
  // plausible-looking but fake result.
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
